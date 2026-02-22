use std::collections::VecDeque;
const INVALID: u8 = 0xFF;

fn base_code(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => INVALID,
    }
}

/// Reverse complement of a 2-bit encoded k-mer of length `n`
#[inline]
fn revcomp(mut kmer: u64, n: usize) -> u64 {
    kmer = ((kmer & 0xCCCC_CCCC_CCCC_CCCC) >> 2) | ((kmer & 0x3333_3333_3333_3333) << 2);
    kmer = ((kmer & 0xF0F0_F0F0_F0F0_F0F0) >> 4) | ((kmer & 0x0F0F_0F0F_0F0F_0F0F) << 4);
    kmer = ((kmer & 0xFF00_FF00_FF00_FF00) >> 8) | ((kmer & 0x00FF_00FF_00FF_00FF) << 8);
    kmer = ((kmer & 0xFFFF_0000_FFFF_0000) >> 16) | ((kmer & 0x0000_FFFF_0000_FFFF) << 16);
    kmer = (kmer >> 32) | (kmer << 32);
    ((!kmer) >> (64 - n * 2)) & ((1u64 << (n * 2)) - 1)
}

#[inline]
fn canonical(kmer: u64, n: usize) -> u64 {
    let rc = revcomp(kmer, n);
    kmer.min(rc)
}

pub struct MinimizerScanner {
    k: usize,
    l: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
    lmer_mask: u64,
}

impl MinimizerScanner {
    pub fn new(k: usize, l: usize, spaced_seed_mask: u64, toggle_mask: u64) -> Self {
        assert!(l <= 31, "l must be <= 31");
        assert!(k >= l, "k must be >= l");
        let lmer_mask = (1u64 << (l * 2)) - 1;
        Self {
            k,
            l,
            spaced_seed_mask,
            toggle_mask: toggle_mask & lmer_mask,
            lmer_mask,
        }
    }

    /// Extract all minimizers from a sequence into a reusable buffer.
    /// Clears `out` before filling it. Returns one minimizer per k-mer window position.
    pub fn scan_into(&self, seq: &[u8], out: &mut Vec<u64>) {
        out.clear();

        if seq.len() < self.k {
            return;
        }

        out.reserve(seq.len() - self.k + 1);
        let mut queue: VecDeque<(u64, usize)> = VecDeque::new();
        let mut lmer: u64 = 0;
        let mut valid: usize = 0;

        for (i, &base) in seq.iter().enumerate() {
            let code = base_code(base);
            if code == INVALID {
                lmer = 0;
                valid = 0;
                queue.clear();
                continue;
            }

            lmer = ((lmer << 2) | code as u64) & self.lmer_mask;
            valid += 1;

            if valid < self.l {
                continue;
            }

            // Compute candidate minimizer value
            let mut can = canonical(lmer, self.l);
            if self.spaced_seed_mask != 0 {
                can &= self.spaced_seed_mask;
            }
            let candidate = can ^ self.toggle_mask;

            // Maintain ascending deque (sliding window minimum)
            while queue.back().map_or(false, |&(c, _)| c > candidate) {
                queue.pop_back();
            }
            let pos = i + 1 - self.l;
            queue.push_back((candidate, pos));

            // Remove expired entries
            let window_start = if pos + self.l >= self.k {
                pos + self.l - self.k
            } else {
                0
            };
            while queue.front().map_or(false, |&(_, p)| p < window_start) {
                queue.pop_front();
            }

            // Emit minimizer if we have a full k-mer window
            if i + 1 >= self.k {
                if let Some(&(min_val, _)) = queue.front() {
                    out.push(min_val ^ self.toggle_mask);
                }
            }
        }
    }

    /// Extract all minimizers from a sequence (allocating version, kept for build path).
    pub fn scan(&self, seq: &[u8]) -> Vec<u64> {
        let mut out = Vec::new();
        self.scan_into(seq, &mut out);
        out
    }

    /// Extract the first minimizer from a sequence (for building).
    pub fn first_minimizer(&self, seq: &[u8]) -> Option<u64> {
        let mins = self.scan(seq);
        mins.into_iter().next()
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn l(&self) -> usize {
        self.l
    }

    pub fn spaced_seed_mask(&self) -> u64 {
        self.spaced_seed_mask
    }

    pub fn toggle_mask(&self) -> u64 {
        self.toggle_mask
    }
}

/// spaced seed mask from a bit pattern string like "1111111111111111111110101010101"
pub fn create_spaced_seed_mask(pattern: &str) -> u64 {
    let mut mask: u64 = 0;
    for (pos, ch) in pattern.chars().rev().enumerate() {
        if ch == '1' {
            mask |= 3u64 << (2 * pos);
        }
    }
    mask
}