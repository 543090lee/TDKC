#[inline(always)]
pub fn extract_accession(id_bytes: &[u8]) -> &str {
    let id_full = std::str::from_utf8(id_bytes).unwrap_or("");
    id_full.split_whitespace().next().unwrap_or("")
}

#[inline]
pub fn segment_ranges(seq_len: usize, k: usize, seg_target_len: usize) -> Vec<(usize, usize)> {
    if seq_len <= seg_target_len {
        return vec![(0, seq_len)];
    }

    let overlap = k.saturating_sub(1);
    let stride = seg_target_len - overlap;
    let mut ranges = Vec::new();
    let mut start = 0;

    while start < seq_len {
        let end = (start + seg_target_len).min(seq_len);
        ranges.push((start, end));
        if end == seq_len {
            break;
        }
        start += stride;
    }

    ranges
}


pub fn init_thread_pool(threads: usize) {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .ok();
}