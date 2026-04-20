// Being used for dequue algorithm for minimizer
pub struct RingDeque {
    vals: [(u64, u32); 64],
    head: u32,
    len: u32,
}

impl RingDeque {
    #[inline(always)]
    pub fn new() -> Self {
        Self {
            vals: [(0u64, 0u32); 64],
            head: 0,
            len: 0,
        }
    }

    #[inline(always)]
    pub fn clear(&mut self) {
        self.head = 0;
        self.len = 0;
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline(always)]
    pub fn back(&self) -> (u64, u32) {
        let idx = (self.head + self.len - 1) & 63;
        self.vals[idx as usize]
    }

    #[inline(always)]
    pub fn front(&self) -> (u64, u32) {
        self.vals[self.head as usize]
    }

    #[inline(always)]
    pub fn pop_back(&mut self) {
        self.len -= 1;
    }

    #[inline(always)]
    pub fn pop_front(&mut self) {
        self.head = (self.head + 1) & 63;
        self.len -= 1;
    }

    #[inline(always)]
    pub fn push_back(&mut self, val: (u64, u32)) {
        let idx = (self.head + self.len) & 63;
        self.vals[idx as usize] = val;
        self.len += 1;
    }
}

#[inline(always)]
pub fn extract_accession(id_bytes: &[u8]) -> &str {
    let id_full = std::str::from_utf8(id_bytes).unwrap_or("");
    id_full.split_whitespace().next().unwrap_or("")
}

pub fn init_thread_pool(threads: usize) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .ok();
}