use std::fs::File;
use std::hash::{BuildHasher, Hasher};
use std::io::{BufWriter, Read};
use std::sync::atomic::{AtomicU8, Ordering};
use std::time::Instant;

use anyhow::{Context, Result};
use hyperloglogplus::{HyperLogLog, HyperLogLogPlus};
use rayon::prelude::*;

use crate::database::DomainBloomFilter;
use crate::minimizer::MinimizerScanner;

#[derive(Default)]
pub struct Fmix64Hasher(u64);

impl Hasher for Fmix64Hasher {
    #[inline(always)]
    fn write_u64(&mut self, i: u64) {
        let mut k = i;
        k ^= k >> 33;
        k = k.wrapping_mul(0xff51afd7ed558ccd);
        k ^= k >> 33;
        k = k.wrapping_mul(0xc4ceb9fe1a85ec53);
        k ^= k >> 33;
        self.0 = k;
    }

    #[inline(always)]
    fn write(&mut self, bytes: &[u8]) {
        let mut arr = [0u8; 8];
        arr.copy_from_slice(&bytes[..8]);
        self.write_u64(u64::from_ne_bytes(arr));
    }

    #[inline(always)]
    fn finish(&self) -> u64 {
        self.0
    }
}

#[derive(Default, Clone)]
pub struct BuildFmix64Hasher;

impl BuildHasher for BuildFmix64Hasher {
    type Hasher = Fmix64Hasher;
    fn build_hasher(&self) -> Self::Hasher {
        Fmix64Hasher::default()
    }
}

struct AtomicBloom {
    data: Vec<AtomicU8>,
    num_hashes: u32,
    num_bits: u64,
}

impl AtomicBloom {
    fn new(expected_items: usize, fpr: f64) -> Self {
        let n = expected_items as f64;
        let ln2 = std::f64::consts::LN_2;
        let m = (-n * fpr.ln() / (ln2 * ln2)).ceil() as u64;
        let k = ((m as f64 / n) * ln2).round() as u32;
        let k = k.max(1);

        let num_bytes = ((m + 7) / 8) as usize;
        let mut data = Vec::with_capacity(num_bytes);
        for _ in 0..num_bytes {
            data.push(AtomicU8::new(0));
        }

        eprintln!(
            "  AtomicBloom: {} bits ({:.2} GiB), {} hashes for {} items @ {:.1}% FPR",
            m,
            num_bytes as f64 / 1_073_741_824.0,
            k,
            expected_items,
            fpr * 100.0
        );

        Self {
            data,
            num_hashes: k,
            num_bits: m,
        }
    }

    #[inline]
    fn insert(&self, item: u64) {
        let (h1, h2) = self.hash_pair(item);
        for i in 0..self.num_hashes {
            let bit_pos = (h1.wrapping_add((i as u64).wrapping_mul(h2))) % self.num_bits;
            let byte_idx = (bit_pos / 8) as usize;
            let bit_mask = 1u8 << (bit_pos % 8);
            self.data[byte_idx].fetch_or(bit_mask, Ordering::Relaxed);
        }
    }

    #[inline(always)]
    fn hash_pair(&self, item: u64) -> (u64, u64) {
        let mut h1 = item;
        h1 ^= h1 >> 33;
        h1 = h1.wrapping_mul(0xff51afd7ed558ccd);
        h1 ^= h1 >> 33;
        h1 = h1.wrapping_mul(0xc4ceb9fe1a85ec53);
        h1 ^= h1 >> 33;

        let mut h2 = item;
        h2 ^= h2 >> 31;
        h2 = h2.wrapping_mul(0x85ebca6b);
        h2 ^= h2 >> 13;
        h2 = h2.wrapping_mul(0xc2b2ae35);
        h2 ^= h2 >> 16;
        h2 |= 1;

        (h1, h2)
    }

    fn into_bytes(self) -> Vec<u8> {
        self.data.into_iter().map(|a| a.into_inner()).collect()
    }

    fn num_hashes(&self) -> u32 {
        self.num_hashes
    }

    fn num_bits(&self) -> u64 {
        self.num_bits
    }
}

pub struct BuildDomainConfig {
    pub db_prefix: String,
    pub threads: usize,
    pub bacteria: Option<String>,
    pub archaea: Option<String>,
    pub viral: Option<String>,
    pub fungi: Option<String>,
}

struct DomainTask {
    name: String,
    fasta_path: String,
    output_ext: String,
}

pub fn run_build_domain(config: BuildDomainConfig) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok();

    let total_start = Instant::now();

    eprintln!("Step 1: Reading Scanner Parameters from Target DB...");
    let (k, l, spaced_seed_mask, toggle_mask) = load_db_meta(&config.db_prefix)?;
    eprintln!(
        "  k={}, l={}, spaced_mask={:016x}, toggle_mask={:016x}",
        k, l, spaced_seed_mask, toggle_mask
    );

    let scanner = MinimizerScanner::new(k, l, spaced_seed_mask, toggle_mask);

    let mut tasks = Vec::new();
    if let Some(path) = config.bacteria { tasks.push(DomainTask { name: "Bacteria".into(), fasta_path: path, output_ext: "bacteria.bloom".into() }); }
    if let Some(path) = config.archaea { tasks.push(DomainTask { name: "Archaea".into(), fasta_path: path, output_ext: "archaea.bloom".into() }); }
    if let Some(path) = config.viral { tasks.push(DomainTask { name: "Viral".into(), fasta_path: path, output_ext: "viral.bloom".into() }); }
    if let Some(path) = config.fungi { tasks.push(DomainTask { name: "Fungi".into(), fasta_path: path, output_ext: "fungi.bloom".into() }); }

    if tasks.is_empty() {
        eprintln!("No domains specified to build. Exiting.");
        return Ok(());
    }

    for task in tasks {
        build_single_domain(&task, &scanner, &config.db_prefix)?;
    }

    eprintln!(
        "\nAll Domain Filters Built! Total time: {:.2}s",
        total_start.elapsed().as_secs_f64()
    );

    Ok(())
}

fn build_single_domain(task: &DomainTask, scanner: &MinimizerScanner, db_prefix: &str) -> Result<()> {
    let start = Instant::now();
    eprintln!("\n========================================");
    eprintln!("Processing Domain: {}", task.name);

    let (unique_count, total_seqs) = estimate_unique_minimizers(&task.fasta_path, scanner)?;

    if unique_count == 0 {
        eprintln!("  WARNING: No valid minimizers found. Skipping.");
        return Ok(());
    }

    eprintln!("\nPass 2: Building Bloom filter with lock-free parallel insertion...");

    let bloom = AtomicBloom::new(unique_count, 0.01);

    let bloom_ref = &bloom;
    let (seq_tx, seq_rx) = crossbeam_channel::bounded::<Vec<Vec<u8>>>(16);
    let fasta_path = task.fasta_path.clone();

    let reader_handle = std::thread::spawn(move || -> Result<()> {
        let mut reader = needletail::parse_fastx_file(&fasta_path)
            .map_err(|e| anyhow::anyhow!("Cannot open FASTA {}: {}", fasta_path, e))?;

        let mut batch = Vec::with_capacity(1000);
        let mut batch_bytes = 0usize;
        let mut current_seq = 0u64;
        let pass2_start = Instant::now();

        while let Some(record) = reader.next() {
            let rec = record.map_err(|e| anyhow::anyhow!("FASTA error: {}", e))?;
            let seq = rec.seq().into_owned();
            batch_bytes += seq.len();
            batch.push(seq);
            current_seq += 1;

            if current_seq % 100_000 == 0 {
                let pct = (current_seq as f64 / total_seqs as f64) * 100.0;
                let elapsed = pass2_start.elapsed().as_secs_f64();
                let rate = current_seq as f64 / elapsed;
                eprint!("\r  [Pass 2] Progress: {} / {} seqs ({:.1}%) | {:.0} seqs/sec   ", current_seq, total_seqs, pct, rate);
            }

            if batch_bytes >= 8 * 1024 * 1024 {
                if seq_tx.send(std::mem::replace(&mut batch, Vec::with_capacity(1000))).is_err() { break; }
                batch_bytes = 0;
            }
        }
        if !batch.is_empty() { let _ = seq_tx.send(batch); }
        eprintln!();
        Ok(())
    });

    let total_inserted = std::sync::atomic::AtomicU64::new(0);
    let total_inserted_ref = &total_inserted;

    seq_rx.into_iter().par_bridge().for_each(|seq_batch| {
        let mut local_mins = Vec::new();
        let mut local_count = 0u64;
        for seq in &seq_batch {
            scanner.scan_into(seq, &mut local_mins);
            for &m in &local_mins {
                if m != u64::MAX {
                    bloom_ref.insert(m);
                    local_count += 1;
                }
            }
        }
        total_inserted_ref.fetch_add(local_count, Ordering::Relaxed);
    });

    match reader_handle.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("Reader thread panicked"),
    }

    let inserted = total_inserted.load(Ordering::Relaxed);
    eprintln!("  Inserted {} total minimizers (including duplicates, which are free for Bloom).", inserted);

    let num_hashes = bloom.num_hashes();
    let num_bits = bloom.num_bits();
    let domain_filter = DomainBloomFilter::new(bloom.into_bytes(), num_hashes, num_bits);

    let popcount = domain_filter.popcount();
    let fill_ratio = popcount as f64 / domain_filter.num_bits as f64;
    let est_fpr = domain_filter.estimated_fpr();
    eprintln!(
        "  Bloom stats: {}/{} bits set ({:.1}% fill), estimated FPR: {:.4}%",
        popcount, domain_filter.num_bits, fill_ratio * 100.0, est_fpr * 100.0
    );

    let out_path = format!("{}.{}", db_prefix, task.output_ext);
    eprintln!("Saving domain filter to {}...", out_path);

    let f = File::create(&out_path)?;
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, f);
    bincode::serialize_into(&mut writer, &domain_filter)?;

    let file_size = std::fs::metadata(&out_path).map(|m| m.len()).unwrap_or(0);
    eprintln!(
        "  Saved {} ({:.2} GiB) in {:.2}s",
        out_path, file_size as f64 / 1_073_741_824.0, start.elapsed().as_secs_f64()
    );

    Ok(())
}

fn estimate_unique_minimizers(fasta_path: &str, scanner: &MinimizerScanner) -> Result<(usize, u64)> {
    eprintln!("Pass 1: Estimating unique minimizers with HyperLogLog++...");
    let start = Instant::now();
    let (seq_tx, seq_rx) = crossbeam_channel::bounded::<Vec<Vec<u8>>>(32);
    let fasta_path_owned = fasta_path.to_string();

    let reader_handle = std::thread::spawn(move || -> Result<u64> {
        let mut reader = needletail::parse_fastx_file(&fasta_path_owned)
            .map_err(|e| anyhow::anyhow!("Cannot open FASTA for HLL {}: {}", fasta_path_owned, e))?;

        let mut batch = Vec::with_capacity(1000);
        let mut batch_bytes = 0;
        let mut seq_count = 0u64;

        while let Some(record) = reader.next() {
            let rec = record.map_err(|e| anyhow::anyhow!("FASTA error: {}", e))?;
            let seq = rec.seq().into_owned();

            batch_bytes += seq.len();
            batch.push(seq);
            seq_count += 1;

            if batch_bytes >= 5 * 1024 * 1024 {
                if seq_tx.send(std::mem::replace(&mut batch, Vec::with_capacity(1000))).is_err() { break; }
                batch_bytes = 0;
            }
        }
        if !batch.is_empty() { let _ = seq_tx.send(batch); }
        Ok(seq_count)
    });

    let mut master_hll = seq_rx.into_iter().par_bridge()
        .fold(
            || HyperLogLogPlus::<u64, BuildFmix64Hasher>::new(16, BuildFmix64Hasher::default()).unwrap(),
            |mut local_hll, seq_batch| {
                let mut local_mins = Vec::new();
                for seq in seq_batch {
                    scanner.scan_into(&seq, &mut local_mins);
                    for &m in &local_mins {
                        if m != u64::MAX { local_hll.insert(&m); }
                    }
                }
                local_hll
            },
        )
        .reduce(
            || HyperLogLogPlus::<u64, BuildFmix64Hasher>::new(16, BuildFmix64Hasher::default()).unwrap(),
            |mut hll_a, hll_b| {
                hll_a.merge(&hll_b).unwrap();
                hll_a
            },
        );

    let seq_count = match reader_handle.join() {
        Ok(Ok(count)) => count,
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("Reader thread panicked"),
    };

    let unique_count = master_hll.count().round() as usize;

    eprintln!("\n  Pass 1 Complete: {} sequences scanned in {:.2}s", seq_count, start.elapsed().as_secs_f64());
    eprintln!("  Estimated Unique Minimizers: {}", unique_count);

    Ok((unique_count, seq_count))
}

fn load_db_meta(prefix: &str) -> Result<(usize, usize, u64, u64)> {
    let mut f = File::open(format!("{}.meta", prefix)).context("Cannot open .meta file to read scanner parameters")?;
    let mut buf4 = [0u8; 4];
    let mut buf8 = [0u8; 8];

    f.read_exact(&mut buf4)?; let k = u32::from_le_bytes(buf4) as usize;
    f.read_exact(&mut buf4)?; let l = u32::from_le_bytes(buf4) as usize;
    f.read_exact(&mut buf8)?; let spaced_seed_mask = u64::from_le_bytes(buf8);
    f.read_exact(&mut buf8)?; let toggle_mask = u64::from_le_bytes(buf8);

    Ok((k, l, spaced_seed_mask, toggle_mask))
}