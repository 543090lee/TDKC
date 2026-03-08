use crate::utils::RingDeque;
const INVALID: u8 = 0xFF;
const BASE_LUT: [u8; 256] = {
    let mut lut = [INVALID; 256];
    lut[b'A' as usize] = 0; lut[b'a' as usize] = 0;
    lut[b'C' as usize] = 1; lut[b'c' as usize] = 1;
    lut[b'G' as usize] = 2; lut[b'g' as usize] = 2;
    lut[b'T' as usize] = 3; lut[b't' as usize] = 3;
    lut
};

#[inline(always)]
fn base_code(c: u8) -> u8 {
    BASE_LUT[c as usize]
}

pub struct MinimizerScanner {
    k: usize,
    l: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
    lmer_mask: u64,
    rc_shift: usize,
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
            rc_shift: (l - 1) * 2,
        }
    }

    // here i am extracting rc the same way as kraken2, in 5 ops
    pub fn scan_into(&self, seq: &[u8], out: &mut Vec<u64>) {
        out.clear();

        if seq.len() < self.k {
            return;
        }

        out.reserve(seq.len() - self.k + 1);
        let mut queue = RingDeque::new();

        let k = self.k;
        let l = self.l;
        let lmer_mask = self.lmer_mask;
        let spaced_seed_mask = self.spaced_seed_mask;
        let toggle_mask = self.toggle_mask;
        let rc_shift = self.rc_shift;
        let has_spaced_seed = spaced_seed_mask != 0;

        let mut lmer: u64 = 0;
        let mut rc_lmer: u64 = 0;
        let mut valid: usize = 0;

        for (i, &base) in seq.iter().enumerate() {
            let code = base_code(base);
            if code == INVALID {
                lmer = 0;
                rc_lmer = 0;
                valid = 0;
                queue.clear();
                continue;
            }

            lmer = ((lmer << 2) | code as u64) & lmer_mask;

            //here complement of code is 3- code
            rc_lmer = (rc_lmer >> 2) | (((3 - code) as u64) << rc_shift);

            valid += 1;

            if valid < l {
                continue;
            }

            let can_lmer = lmer.min(rc_lmer);

            let mut can = can_lmer;
            if has_spaced_seed {
                can &= spaced_seed_mask;
            }
            let candidate = can ^ toggle_mask;

            // I am using deque algorithm here. the queue is always in the descending order.
            // when you are adding a candidate to the queue, compare it to the back, and if it's smaller
            // then we know that the back will never be the window minimum, pop it. then check the next back again
            // if finally bigger, then you just add it.
            // so it's amortized O(1), 
            while !queue.is_empty() && queue.back().0 > candidate {
                queue.pop_back();
            }
            let pos = (i + 1 - l) as u32;
            queue.push_back((candidate, pos));

            let window_start = if (pos as usize) + l >= k {
                (pos as usize + l - k) as u32
            } else {
                0
            };

            // remove the front, if it now falls outside the window. this might be the minimizer too
            while !queue.is_empty() && queue.front().1 < window_start {
                queue.pop_front();
            }

            // Emit minimizer when we have a full k-mer window
            if i + 1 >= k {
                let (min_val, _) = queue.front();
                out.push(min_val ^ toggle_mask);
            }
        }
    }

    pub fn k(&self) -> usize {
        self.k
    }
    
}

// spaced seed mask from a bit pattern string like "1111111111111111111110101010101"
pub fn create_spaced_seed_mask(pattern: &str) -> u64 {
    let mut mask: u64 = 0;
    for (pos, ch) in pattern.chars().rev().enumerate() {
        if ch == '1' {
            mask |= 3u64 << (2 * pos);
        }
    }
    mask
}
