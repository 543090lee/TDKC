
pub fn compression_stats(&self) -> CompressionStats {

    let mut class_sizes: Vec<u32> = Vec::with_capacity(self.num_classes.saturating_sub(1));
    let mut total_acc_refs: u64 = 0;

    for cid in 1..self.num_classes {
        let mut pos = self.offsets[cid] as usize;
        let count = vbyte_decode(&self.blob, &mut pos) as u32;
        class_sizes.push(count);
        total_acc_refs += count as u64;
    }

    let num_real_classes = class_sizes.len();
    class_sizes.sort_unstable();

    let (min_size, max_size, median_size, mean_size) = if num_real_classes == 0 {
        (0u32, 0u32, 0.0f64, 0.0f64)
    } else {
        let min_s = class_sizes[0];
        let max_s = *class_sizes.last().unwrap();
        let median = if num_real_classes % 2 == 0 {
            (class_sizes[num_real_classes / 2 - 1] as f64
                + class_sizes[num_real_classes / 2] as f64) / 2.0
        } else {
            class_sizes[num_real_classes / 2] as f64
        };
        let mean = total_acc_refs as f64 / num_real_classes as f64;
        (min_s, max_s, median, mean)
    };

    let pct = |p: f64| -> u32 {
        if num_real_classes == 0 { return 0; }
        let idx = ((num_real_classes as f64 - 1.0) * p).round() as usize;
        class_sizes[idx.min(num_real_classes - 1)]
    };
    let p25 = pct(0.25);
    let p75 = pct(0.75);
    let p90 = pct(0.90);
    let p99 = pct(0.99);

    let mut hist = [0usize; 8];
    for &s in &class_sizes {
        let bucket = match s {
            0 => continue,
            1 => 0,
            2 => 1,
            3..=5 => 2,
            6..=10 => 3,
            11..=50 => 4,
            51..=100 => 5,
            101..=1000 => 6,
            _ => 7,
        };
        hist[bucket] += 1;
    }

    let singletons = hist[0];
    let shared_classes = num_real_classes - singletons;

    let mut minimizers_pointing_to_singleton: u64 = 0;
    let mut minimizers_pointing_to_shared: u64 = 0;

    let mut size_by_cid: Vec<u32> = vec![0u32; self.num_classes];
    for cid in 1..self.num_classes {
        let mut pos = self.offsets[cid] as usize;
        size_by_cid[cid] = vbyte_decode(&self.blob, &mut pos);
    }
    for &cid in &self.dense_ids {
        let sz = size_by_cid[cid as usize];
        if sz == 1 {
            minimizers_pointing_to_singleton += 1;
        } else if sz > 1 {
            minimizers_pointing_to_shared += 1;
        }
    }
    let total_tracked_minimizers = self.dense_ids.len() as u64;

    let mut num_deltas: u64 = 0;
    let mut deltas_1byte: u64 = 0;
    let mut deltas_2byte: u64 = 0;
    let mut deltas_3byte: u64 = 0;
    let mut deltas_4plus: u64 = 0;
    let mut delta_sum: u64 = 0;
    let mut delta_max: u32 = 0;

    for cid in 1..self.num_classes {
        let mut pos = self.offsets[cid] as usize;
        let count = vbyte_decode(&self.blob, &mut pos);
        let mut prev: u32 = 0;
        for _ in 0..count {
            let delta = vbyte_decode(&self.blob, &mut pos);
            prev += delta;
            num_deltas += 1;
            delta_sum += delta as u64;
            if delta > delta_max { delta_max = delta; }
            if delta < 128          { deltas_1byte += 1; }
            else if delta < 16_384  { deltas_2byte += 1; }
            else if delta < 2_097_152 { deltas_3byte += 1; }
            else                    { deltas_4plus += 1; }
        }
    }
    let mean_delta = if num_deltas > 0 { delta_sum as f64 / num_deltas as f64 } else { 0.0 };

    let bitset_bytes = self.bitset.len() * 8;
    let rank_cache_bytes = self.rank_cache.len() * 4;
    let dense_id_width = if self.num_classes <= 256 { 1 }
                        else if self.num_classes <= 65_536 { 2 }
                        else { 4 };
    let dense_ids_bytes = self.dense_ids.len() * dense_id_width;
    let offsets_bytes = self.offsets.len() * 4;
    let blob_bytes = self.blob.len();
    let total_actual_bytes = bitset_bytes + rank_cache_bytes + dense_ids_bytes + offsets_bytes + blob_bytes;

    let mut naive_per_minimizer_bytes: u64 = 0;
    for &cid in &self.dense_ids {
        let sz = size_by_cid[cid as usize] as u64;
        naive_per_minimizer_bytes += 4 + 4 * sz;
    }
    
    let naive_index_bytes = (self.num_slots * 4) as u64;
    let naive_total_bytes = naive_per_minimizer_bytes + naive_index_bytes;
    let mut eqclass_only_blob_bytes: u64 = 0;
    for cid in 1..self.num_classes {
        let sz = size_by_cid[cid as usize] as u64;
        eqclass_only_blob_bytes += 4 + 4 * sz;
    }
    let eqclass_only_total_bytes = bitset_bytes as u64
        + rank_cache_bytes as u64
        + dense_ids_bytes as u64
        + offsets_bytes as u64
        + eqclass_only_blob_bytes;

    CompressionStats {
        num_slots: self.num_slots,
        total_tracked_minimizers,
        num_classes: num_real_classes,
        compression_ratio: if num_real_classes > 0 {
            total_tracked_minimizers as f64 / num_real_classes as f64
        } else { 0.0 },
        total_accession_refs: total_acc_refs,
        min_size, max_size, median_size, mean_size,
        p25, p75, p90, p99,
        hist,
        singletons,
        shared_classes,
        minimizers_pointing_to_singleton,
        minimizers_pointing_to_shared,
        num_deltas,
        mean_delta,
        delta_max,
        deltas_1byte, deltas_2byte, deltas_3byte, deltas_4plus,
        bitset_bytes,
        rank_cache_bytes,
        dense_ids_bytes,
        dense_id_width,
        offsets_bytes,
        blob_bytes,
        total_actual_bytes,
        naive_per_minimizer_bytes,
        naive_index_bytes,
        naive_total_bytes,
        eqclass_only_blob_bytes,
        eqclass_only_total_bytes,
    }
}