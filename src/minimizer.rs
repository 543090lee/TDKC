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

pub const SPACED_PATTERN: &str = "1111111111111111101010101010101";
pub const TOGGLE_MASK: u64 = 0xe37e28c4271b5a2d;

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
        assert!(l <= 31, "l too big! l has to be <= 31");
        assert!(k >= l, "k too small! k has to be >= l");
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

    pub fn scan_into(&self, seq: &[u8], out: &mut Vec<u64>) {
        out.clear();

        if seq.len() < self.k {
            return;
        }

        let num_kmers = seq.len() - self.k + 1;
        out.reserve(num_kmers);
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

                if i + 1 >= k {
                    out.push(u64::MAX);
                }
                continue;
            }

            lmer = ((lmer << 2) | code as u64) & lmer_mask;
            rc_lmer = (rc_lmer >> 2) | (((3 - code) as u64) << rc_shift);
            valid += 1;

            if valid >= l {
                let can_lmer = lmer.min(rc_lmer);
                let mut can = can_lmer;
                if has_spaced_seed {
                    can &= spaced_seed_mask;
                }
                let candidate = can ^ toggle_mask;

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

                while !queue.is_empty() && queue.front().1 < window_start {
                    queue.pop_front();
                }
            }

            if i + 1 < k {
                continue;
            }

            if valid >= k - 1 && !queue.is_empty() {
                let (min_val, _) = queue.front();
                out.push(min_val ^ toggle_mask);
            } else {
                out.push(u64::MAX);
            }
        }

        debug_assert_eq!(
            out.len(),
            num_kmers,
            "scan_into output length {} != expected {} for seq len {}",
            out.len(),
            num_kmers,
            seq.len()
        );
    }

    pub fn k(&self) -> usize {
        self.k
    }
}

// It has the same configuration as Kraken2's default space seed mask pattern
pub fn create_spaced_seed_mask(pattern: &str) -> u64 {
    let mut mask: u64 = 0;
    for (pos, ch) in pattern.chars().rev().enumerate() {
        if ch == '1' {
            mask |= 3u64 << (2 * pos);
        }
    }
    mask
}

struct RingDeque {
    vals: [(u64, u32); 64],
    head: u32,
    len: u32,
}

impl RingDeque {
    #[inline(always)]
    fn new() -> Self {
        Self {
            vals: [(0u64, 0u32); 64],
            head: 0,
            len: 0,
        }
    }

    #[inline(always)]
    fn clear(&mut self) {
        self.head = 0;
        self.len = 0;
    }

    #[inline(always)]
    fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline(always)]
    fn back(&self) -> (u64, u32) {
        let idx = (self.head + self.len - 1) & 63;
        self.vals[idx as usize]
    }

    #[inline(always)]
    fn front(&self) -> (u64, u32) {
        self.vals[self.head as usize]
    }

    #[inline(always)]
    fn pop_back(&mut self) {
        self.len -= 1;
    }

    #[inline(always)]
    fn pop_front(&mut self) {
        self.head = (self.head + 1) & 63;
        self.len -= 1;
    }

    #[inline(always)]
    fn push_back(&mut self, val: (u64, u32)) {
        let idx = (self.head + self.len) & 63;
        self.vals[idx as usize] = val;
        self.len += 1;
    }
}
