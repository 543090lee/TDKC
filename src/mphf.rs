
use anyhow::{bail, Context, Result};
use rayon::prelude::*;

#[repr(C, align(64))]
#[derive(Clone, Copy)]
pub struct Block(pub [u64; 8]);

pub const WORDS_PER_BLOCK: usize = 7;

pub struct Level {
    nbits: u64,
    words: Vec<u64>,
    ranks: Vec<u64>,
}

pub struct FastMphf {
    levels: Vec<Level>,
    total_keys: u64,
}

const P0: u64 = 0xa076_1d64_78bd_642f;
const P1: u64 = 0xe703_7ed1_a0b4_28db;
const P5: u64 = 0xeb44_acca_b455_d165;

#[inline(always)]
fn wymum(a: u64, b: u64) -> u64 {
    let r = (a as u128) * (b as u128);
    ((r >> 64) ^ r) as u64
}

#[inline(always)]
pub fn boomphf_hash(level: u64, key: u64) -> u64 {
    let seed = 1u64.wrapping_shl((level + level) as u32);
    let swapped = (key << 32) | (key >> 32);
    let h = wymum(seed ^ P0, swapped ^ P1);
    wymum(h, 8 ^ P5)
}

#[inline(always)]
fn fold(v: u64) -> u32 {
    ((v & 0xFFFF_FFFF) as u32) ^ ((v >> 32) as u32)
}

#[inline(always)]
pub fn hashmod(level: u64, key: u64, nbits: u64) -> u64 {
    let h = boomphf_hash(level, key);
    if nbits < (1u64 << 32) {
        ((fold(h) as u64) * nbits) >> 32
    } else {
        h % nbits
    }
}

fn read_u64s<R: std::io::Read>(r: &mut R, n: usize, what: &str) -> Result<Vec<u64>> {
    let mut v = vec![0u64; n];
    let dst = unsafe { std::slice::from_raw_parts_mut(v.as_mut_ptr() as *mut u8, n * 8) };
    r.read_exact(dst).with_context(|| format!("truncated .mphf file ({})", what))?;
    #[cfg(target_endian = "big")]
    for x in v.iter_mut() {
        *x = u64::from_le(*x);
    }
    Ok(v)
}

impl FastMphf {
    #[allow(dead_code)]
    pub fn from_bincode_bytes(bytes: &[u8]) -> Result<Self> {
        let mut cur = std::io::Cursor::new(bytes);
        Self::from_reader(&mut cur, bytes.len() as u64)
    }

    pub fn from_reader<R: std::io::Read>(r: &mut R, total_len: u64) -> Result<Self> {
        let mut consumed: u64 = 0;
        let rd8 = |r: &mut R, consumed: &mut u64| -> Result<u64> {
            let mut b = [0u8; 8];
            r.read_exact(&mut b).context("truncated .mphf file")?;
            *consumed += 8;
            Ok(u64::from_le_bytes(b))
        };
        let nlevels = rd8(r, &mut consumed)? as usize;
        if nlevels == 0 || nlevels > 256 {
            bail!("unexpected number of MPHF levels: {}", nlevels);
        }
        let mut levels = Vec::with_capacity(nlevels);
        let mut total_keys: u64 = 0;
        for li in 0..nlevels {
            let nbits = rd8(r, &mut consumed)?;
            let nwords = rd8(r, &mut consumed)? as usize;
            if nwords != ((nbits + 63) / 64) as usize {
                bail!("MPHF level {}: {} words for {} bits", li, nwords, nbits);
            }
            let words = read_u64s(r, nwords, "words")?;
            consumed += (nwords * 8) as u64;
            let nranks = rd8(r, &mut consumed)? as usize;
            if nranks != (nwords + 7) / 8 {
                bail!("MPHF level {}: {} rank entries for {} words", li, nranks, nwords);
            }
            let ranks = read_u64s(r, nranks, "ranks")?;
            consumed += (nranks * 8) as u64;
            total_keys += words.par_iter().map(|w| w.count_ones() as u64).sum::<u64>();
            levels.push(Level { nbits, words, ranks });
        }
        if consumed != total_len {
            bail!(".mphf has {} trailing bytes after {} levels", total_len as i64 - consumed as i64, nlevels);
        }
        Ok(Self { levels, total_keys })
    }

    pub fn load(path: &str) -> Result<Self> {
        let f = std::fs::File::open(path).with_context(|| format!("Cannot open {}", path))?;
        let len = f.metadata()?.len();
        let mut br = std::io::BufReader::with_capacity(8 << 20, f);
        let r = Self::from_reader(&mut br, len);
        #[cfg(target_os = "linux")]
        if let Ok(ref m) = r {
            for lvl in &m.levels {
                unsafe {
                    libc::madvise(lvl.words.as_ptr() as *mut libc::c_void, lvl.words.len() * 8, libc::MADV_HUGEPAGE);
                }
            }
        }
        r
    }

    #[inline(always)]
    pub fn num_levels(&self) -> usize {
        self.levels.len()
    }

    pub fn total_keys(&self) -> u64 {
        self.total_keys
    }

    #[inline(always)]
    pub fn level_index(&self, level: usize, key: u64) -> u64 {
        hashmod(level as u64, key, self.levels[level].nbits)
    }

    #[inline(always)]
    pub fn prefetch(&self, level: usize, idx: u64) {
        let word = (idx >> 6) as usize;
        let lvl = &self.levels[level];
        unsafe {
            prefetch_read(lvl.words.as_ptr().add(word));
            prefetch_read(lvl.ranks.as_ptr().add(word >> 3));
        }
    }

    #[inline(always)]
    pub fn test_and_rank(&self, level: usize, idx: u64) -> Option<u64> {
        let word = (idx >> 6) as usize;
        let bit = (idx & 63) as u32;
        let lvl = &self.levels[level];
        let w = lvl.words[word];
        if (w >> bit) & 1 == 0 {
            return None;
        }
        let mut r = lvl.ranks[word >> 3];
        for t in (word & !7)..word {
            r += lvl.words[t].count_ones() as u64;
        }
        if bit > 0 {
            r += (w & ((1u64 << bit) - 1)).count_ones() as u64;
        }
        Some(r)
    }

    #[allow(dead_code)]
    #[inline]
    pub fn try_hash(&self, key: u64) -> Option<u64> {
        for level in 0..self.levels.len() {
            let idx = self.level_index(level, key);
            if let Some(r) = self.test_and_rank(level, idx) {
                return Some(r);
            }
        }
        None
    }
}

#[inline(always)]
pub fn prefetch_read<T>(p: *const T) {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        core::arch::x86_64::_mm_prefetch(p as *const i8, core::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(target_arch = "aarch64")]
    unsafe {
        core::arch::asm!("prfm pldl1keep, [{0}]", in(reg) p, options(nostack, preserves_flags));
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        let _ = p;
    }
}
