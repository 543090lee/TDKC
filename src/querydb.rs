
use crate::database::TaxIdx;
use crate::database::{read_meta, read_names, read_taxmap, DomainBloomFilter, EqClassAccessions, MetaHeader};
use crate::hash::compute_fingerprint;
use crate::mphf::{prefetch_read, Block, FastMphf, WORDS_PER_BLOCK};
use crate::taxonomy::{resolve_domain_lca, TAXID_ARCHAEA, TAXID_BACTERIA, TAXID_CELLULAR, TAXID_FUNGI, TAXID_ROOT, TAXID_VIRAL};
use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use std::fs::File;
use std::path::Path;

pub struct AccIndex {
    blocks: Vec<Block>,
    num_words: usize,
    dense_ids: Vec<u32>,
    offsets: Vec<u32>,
    blob: Vec<u8>,
    num_classes: usize,
}

impl AccIndex {
    fn from_eq(acc: EqClassAccessions) -> Self {
        let (bitset, dense_ids, offsets, blob, num_classes) = acc.into_parts();
        let num_words = bitset.len();
        let nblocks = (num_words + WORDS_PER_BLOCK - 1) / WORDS_PER_BLOCK;
        let mut blocks = vec![Block([0u64; 8]); nblocks.max(1)];
        let mut run = 0u64;
        for (b, blk) in blocks.iter_mut().enumerate() {
            blk.0[7] = run;
            for s in 0..WORDS_PER_BLOCK {
                let w = b * WORDS_PER_BLOCK + s;
                if w < num_words {
                    blk.0[s] = bitset[w];
                    run += bitset[w].count_ones() as u64;
                }
            }
        }
        Self { blocks, num_words, dense_ids, offsets, blob, num_classes }
    }

    #[inline(always)]
    pub fn block_ptr(&self, idx: usize) -> *const Block {
        let b = (idx >> 6) / WORDS_PER_BLOCK;
        unsafe { self.blocks.as_ptr().add(b.min(self.blocks.len() - 1)) }
    }

    #[inline(always)]
    pub fn dense_pos(&self, idx: usize) -> Option<u32> {
        let word = idx >> 6;
        if word >= self.num_words {
            return None;
        }
        let bit = (idx & 63) as u32;
        let b = word / WORDS_PER_BLOCK;
        let s = word - b * WORDS_PER_BLOCK;
        let blk = &self.blocks[b].0;
        let w = blk[s];
        if (w >> bit) & 1 == 0 {
            return None;
        }
        let mut r = blk[7];
        for t in 0..s {
            r += blk[t].count_ones() as u64;
        }
        if bit > 0 {
            r += (w & ((1u64 << bit) - 1)).count_ones() as u64;
        }
        Some(r as u32)
    }

    #[inline(always)]
    pub fn dense_ptr(&self, dense_pos: u32) -> *const u32 {
        unsafe { self.dense_ids.as_ptr().add((dense_pos as usize).min(self.dense_ids.len().saturating_sub(1))) }
    }

    #[inline(always)]
    pub fn class_id(&self, dense_pos: u32) -> u32 {
        self.dense_ids[dense_pos as usize]
    }

    pub fn num_dense(&self) -> usize {
        self.dense_ids.len()
    }

    pub fn decode_class_into(&self, class_id: u32, out: &mut Vec<u32>) {
        out.clear();
        let cid = class_id as usize;
        if cid == 0 || cid >= self.num_classes {
            return;
        }
        let mut pos = self.offsets[cid] as usize;
        let count = crate::compression::vbyte_decode(&self.blob, &mut pos) as usize;
        let mut prev = 0u32;
        for _ in 0..count {
            let delta = crate::compression::vbyte_decode(&self.blob, &mut pos);
            prev += delta;
            out.push(prev);
        }
    }
}

#[derive(Clone, Copy)]
pub struct Lookup {
    pub key: u64,
    h: u64,
    idx: u64,
    level: u16,
    state: u8,
    pub is_hit: bool,
    pub taxid_idx: TaxIdx,
    pub class_id: u32,
    pub dense_pos: u32,
    pub bg_taxid: u32,
}

const ST_PROBING: u8 = 0;
const ST_HAVE_IDX: u8 = 1;
const ST_MISS: u8 = 2;
const ST_HIT: u8 = 3;

impl Lookup {
    #[inline(always)]
    pub fn new(key: u64) -> Self {
        Self { key, h: 0, idx: 0, level: 0, state: ST_PROBING, is_hit: false, taxid_idx: 0, class_id: u32::MAX, dense_pos: 0, bg_taxid: 0 }
    }
}

pub struct QueryDb {
    pub k: usize,
    pub l: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    pub num_minimizers: usize,
    mphf: FastMphf,
    fptax: Vec<u8>,
    stride: usize,
    pub index_to_taxid: Vec<u32>,
    pub ancestor_matrix: Vec<u8>,
    pub accessions: Option<AccIndex>,
    pub names_map: Option<rustc_hash::FxHashMap<u32, String>>,

    pub mask_words: usize,
    pub desc_lists: Vec<Vec<TaxIdx>>,
    pub anc_masks: Vec<u64>,
    pub domain_mask: Vec<u64>,
    pub domain_nodes: Vec<(u32, TaxIdx)>,
}

impl QueryDb {
    pub fn load(prefix: &str, load_accessions: bool) -> Result<Self> {
        eprintln!("Loading database from {}", prefix);
        let t0 = std::time::Instant::now();

        let MetaHeader { k, l, spaced_seed_mask, toggle_mask, num_minimizers, has_acc, taxid_width: meta_width } = read_meta(prefix)?;

        let mphf = FastMphf::load(&format!("{}.mphf", prefix)).context("Cannot load .mphf file")?;
        if mphf.total_keys() != num_minimizers as u64 {
            eprintln!(
                "Warning: .mphf holds {} keys but .meta says {} minimizers",
                mphf.total_keys(),
                num_minimizers
            );
        }

        let fptax = {
            let ffp = File::open(format!("{}.fp", prefix)).context("Cannot open .fp file")?;
            let ftx = File::open(format!("{}.taxid", prefix)).context("Cannot open .taxid file")?;
            if (ffp.metadata()?.len() as usize) < num_minimizers * 2 {
                bail!(".fp file is shorter than expected");
            }
            let tw = crate::database::infer_taxid_width(ftx.metadata()?.len(), num_minimizers)?;
            crate::database::check_taxid_width(meta_width, tw)?;
            let stride = 2 + tw;
            let mut out: Vec<u8> = Vec::with_capacity(num_minimizers * stride + 8);
            unsafe { out.set_len(num_minimizers * stride) };
            const CH: usize = 1 << 20;
            let err = std::sync::Mutex::new(None::<std::io::Error>);
            out.par_chunks_mut(CH * stride).enumerate().for_each(|(ci, dst)| {
                let base = ci * CH;
                let n = dst.len() / stride;
                let mut fp = vec![0u8; n * 2];
                let mut tx = vec![0u8; n * tw];
                let r1 = crate::pgz::read_at(&ffp, &mut fp, (base * 2) as u64);
                let r2 = crate::pgz::read_at(&ftx, &mut tx, (base * tw) as u64);
                if let Err(e) = r1.and(r2) {
                    *err.lock().unwrap() = Some(e);
                    return;
                }
                for i in 0..n {
                    let o = i * stride;
                    dst[o] = fp[i * 2];
                    dst[o + 1] = fp[i * 2 + 1];
                    dst[o + 2] = tx[i * tw];
                    if tw == 2 {
                        dst[o + 3] = tx[i * 2 + 1];
                    }
                }
            });
            if let Some(e) = err.into_inner().unwrap() {
                return Err(e).context("reading .fp/.taxid");
            }
            out.extend_from_slice(&[0u8; 8]);
            #[cfg(target_os = "linux")]
            unsafe {
                libc::madvise(out.as_mut_ptr() as *mut libc::c_void, out.len(), libc::MADV_HUGEPAGE);
            }
            (out, stride)
        };
        let (fptax, stride) = fptax;

        let (index_to_taxid, ancestor_matrix) = read_taxmap(prefix)?;
        if index_to_taxid.len() > crate::database::MAX_TAXA {
            bail!("database has {} target taxa; the taxid index supports at most {}", index_to_taxid.len(), crate::database::MAX_TAXA);
        }
        if stride == 3 && index_to_taxid.len() > 256 {
            bail!(".taxid is 1 byte per minimizer but .taxmap lists {} taxa", index_to_taxid.len());
        }

        let names_map = read_names(prefix);

        let acc_path = format!("{}.accession", prefix);
        let accessions = if load_accessions && has_acc && Path::new(&acc_path).exists() {
            Some(AccIndex::from_eq(EqClassAccessions::load(&acc_path)?))
        } else {
            None
        };

        let n = index_to_taxid.len();
        let is_anc = |a: usize, d: usize| ancestor_matrix[a * n + d] != 0;
        let mask_words = (n + 63) / 64;
        let mut desc_lists = vec![Vec::new(); n];
        let mut anc_masks = vec![0u64; n * mask_words];
        for a in 0..n {
            for d in 0..n {
                if is_anc(a, d) {
                    desc_lists[a].push(d as TaxIdx);
                    anc_masks[d * mask_words + a / 64] |= 1u64 << (a % 64);
                }
            }
        }
        let mut domain_mask = vec![0u64; mask_words];
        let mut domain_nodes = Vec::new();
        for (i, &t) in index_to_taxid.iter().enumerate() {
            if matches!(t, TAXID_ROOT | TAXID_BACTERIA | TAXID_ARCHAEA | TAXID_FUNGI | TAXID_VIRAL | TAXID_CELLULAR) {
                domain_mask[i / 64] |= 1u64 << (i % 64);
            }
        }
        for t in [TAXID_ROOT, TAXID_BACTERIA, TAXID_ARCHAEA, TAXID_FUNGI, TAXID_VIRAL, TAXID_CELLULAR] {
            if let Some(p) = index_to_taxid.iter().position(|&x| x == t) {
                domain_nodes.push((t, p as TaxIdx));
            }
        }

        eprintln!(
            "Loaded {} minimizers, {} target taxa, {} MPHF levels in {:.2}s",
            num_minimizers,
            n,
            mphf.num_levels(),
            t0.elapsed().as_secs_f64()
        );

        Ok(Self {
            k,
            l,
            spaced_seed_mask,
            toggle_mask,
            num_minimizers,
            mphf,
            fptax,
            stride,
            index_to_taxid,
            ancestor_matrix,
            accessions,
            names_map,
            mask_words,
            desc_lists,
            anc_masks,
            domain_mask,
            domain_nodes,
        })
    }

    #[inline(always)]
    pub fn true_taxid(&self, idx: TaxIdx) -> u32 {
        self.index_to_taxid.get(idx as usize).copied().unwrap_or(0)
    }

    #[inline(always)]
    pub fn is_ancestor(&self, ancestor_idx: TaxIdx, descendant_idx: TaxIdx) -> bool {
        let n = self.index_to_taxid.len();
        let a = ancestor_idx as usize;
        let d = descendant_idx as usize;
        if a >= n || d >= n {
            return false;
        }
        self.ancestor_matrix[a * n + d] != 0
    }

    #[inline(always)]
    pub fn get_taxid_index(&self, taxid: u32) -> Option<TaxIdx> {
        for &(t, i) in &self.domain_nodes {
            if t == taxid {
                return Some(i);
            }
        }
        self.index_to_taxid.iter().position(|&t| t == taxid).map(|p| p as TaxIdx)
    }

    pub fn get_taxid_name(&self, taxid: u32) -> String {
        if taxid == u32::MAX {
            return "\t".to_string();
        }
        if let Some(names_map) = &self.names_map {
            if let Some(name) = names_map.get(&taxid) {
                return name.clone();
            }
        }
        "Unknown".to_string()
    }

    #[inline(always)]
    fn fptax_ptr(&self, idx: u64) -> *const u8 {
        unsafe { self.fptax.as_ptr().add((idx as usize) * self.stride) }
    }

    #[inline(always)]
    pub fn anc_mask(&self, i: usize) -> &[u64] {
        &self.anc_masks[i * self.mask_words..(i + 1) * self.mask_words]
    }

    pub fn lookup_batch(&self, slots: &mut [Lookup], scratch: &mut Vec<u32>, bg_filters: &[(DomainBloomFilter, u32)]) {
        let nlevels = self.mphf.num_levels();
        let acc = self.accessions.as_ref();

        scratch.clear();
        for (i, s) in slots.iter_mut().enumerate() {
            s.level = 0;
            s.state = ST_PROBING;
            s.is_hit = false;
            s.taxid_idx = 0;
            s.class_id = u32::MAX;
            s.dense_pos = 0;
            s.bg_taxid = 0;
            s.h = self.mphf.level_index(0, s.key);
            self.mphf.prefetch(0, s.h);
            scratch.push(i as u32);
        }

        let mut level = 0usize;
        while !scratch.is_empty() {
            let mut wr = 0usize;
            for r in 0..scratch.len() {
                let i = scratch[r] as usize;
                let s = &mut slots[i];
                match self.mphf.test_and_rank(level, s.h) {
                    Some(rank) => {
                        s.idx = rank;
                        s.state = ST_HAVE_IDX;
                        if (rank as usize) < self.num_minimizers {
                            prefetch_read(self.fptax_ptr(rank));
                            if let Some(a) = acc {
                                prefetch_read(a.block_ptr(rank as usize));
                            }
                        }
                    }
                    None => {
                        if level + 1 < nlevels {
                            s.level = (level + 1) as u16;
                            s.h = self.mphf.level_index(level + 1, s.key);
                            self.mphf.prefetch(level + 1, s.h);
                            scratch[wr] = i as u32;
                            wr += 1;
                        } else {
                            s.state = ST_MISS;
                        }
                    }
                }
            }
            scratch.truncate(wr);
            level += 1;
        }

        let mut any_acc_hit = false;
        for s in slots.iter_mut() {
            if s.state == ST_HAVE_IDX {
                let idx = s.idx as usize;
                if idx < self.num_minimizers {
                    let p = idx * self.stride;
                    let fp = u16::from_le_bytes([self.fptax[p], self.fptax[p + 1]]);
                    if fp == compute_fingerprint(s.key) {
                        s.state = ST_HIT;
                        s.is_hit = true;
                        s.taxid_idx = if self.stride == 3 {
                            self.fptax[p + 2] as TaxIdx
                        } else {
                            u16::from_le_bytes([self.fptax[p + 2], self.fptax[p + 3]])
                        };
                        if let Some(a) = acc {
                            if let Some(dp) = a.dense_pos(idx) {
                                s.dense_pos = dp;
                                s.class_id = u32::MAX - 1;
                                prefetch_read(a.dense_ptr(dp));
                                any_acc_hit = true;
                            }
                        }
                        continue;
                    }
                }
                s.state = ST_MISS;
            }
        }
        if any_acc_hit {
            let a = acc.unwrap();
            for s in slots.iter_mut() {
                if s.state == ST_HIT && s.class_id == u32::MAX - 1 {
                    s.class_id = a.class_id(s.dense_pos);
                }
            }
        }

        if !bg_filters.is_empty() {
            for s in slots.iter_mut() {
                if s.state == ST_MISS {
                    let mut hit_taxids = [0u32; 4];
                    let mut match_count = 0usize;
                    for (filter, t_id) in bg_filters {
                        if filter.contains(s.key) {
                            if match_count < 4 {
                                hit_taxids[match_count] = *t_id;
                            }
                            match_count += 1;
                        }
                    }
                    s.bg_taxid = match match_count {
                        0 => 0,
                        1 => hit_taxids[0],
                        _ => resolve_domain_lca(hit_taxids[..match_count.min(4)].iter()),
                    };
                }
            }
        }
    }

}
