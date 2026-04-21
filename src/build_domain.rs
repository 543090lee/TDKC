use std::fs::File;
use std::io::{BufWriter, Read};
use std::sync::atomic::{AtomicU8, Ordering};
use std::time::Instant;
use anyhow::{Context, Result};
use hyperloglogplus::{HyperLogLog, HyperLogLogPlus};
use rayon::prelude::*;
use crate::database::DomainBloomFilter;
use crate::minimizer::MinimizerScanner;
use crate::utils::init_thread_pool;
use crate::hash::{BuildFmix64Hasher,bloom_hash_pair};

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

        eprintln!("{} bits ({:.2} GiB), {} hashes for {} items @ {:.1}% FPR", m, num_bytes as f64 / 1_073_741_824.0, k, expected_items, fpr * 100.0);

        Self {
            data,
            num_hashes: k,
            num_bits: m,
        }
    }

    #[inline]
    fn insert(&self, item: u64) {
        let (h1, h2) = bloom_hash_pair(item);
        for i in 0..self.num_hashes {
            let bit_pos = (h1.wrapping_add((i as u64).wrapping_mul(h2))) % self.num_bits;
            let byte_idx = (bit_pos / 8) as usize;
            let bit_mask = 1u8 << (bit_pos % 8);
            self.data[byte_idx].fetch_or(bit_mask, Ordering::Relaxed);
        }
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
    pub fpr: f64,
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
    init_thread_pool(config.threads);

    let total_start = Instant::now();
    let (k, l, spaced_seed_mask, toggle_mask) = load_db_meta(&config.db_prefix)?;
    let scanner = MinimizerScanner::new(k, l, spaced_seed_mask, toggle_mask);

    let mut tasks = Vec::new();
    if let Some(path) = config.bacteria { tasks.push(DomainTask { name: "Bacteria".into(), fasta_path: path, output_ext: "bacteria.bloom".into() }); }
    if let Some(path) = config.archaea { tasks.push(DomainTask { name: "Archaea".into(), fasta_path: path, output_ext: "archaea.bloom".into() }); }
    if let Some(path) = config.viral { tasks.push(DomainTask { name: "Viral".into(), fasta_path: path, output_ext: "viral.bloom".into() }); }
    if let Some(path) = config.fungi { tasks.push(DomainTask { name: "Fungi".into(), fasta_path: path, output_ext: "fungi.bloom".into() }); }

    if tasks.is_empty() {
        eprintln!("No domains specified.");
        return Ok(());
    }

    for task in tasks {
        build_single_domain(&task, &scanner, &config.db_prefix, config.fpr)?;
    }

    eprintln!(
        "\nDone. It took {:.2}s",
        total_start.elapsed().as_secs_f64()
    );
    Ok(())
}

fn build_single_domain(task: &DomainTask, scanner: &MinimizerScanner, db_prefix: &str, fpr: f64) -> Result<()> {
    eprintln!("Processing Domain: {}", task.name);
    let (unique_count, total_seqs) = estimate_unique_minimizers(&task.fasta_path, scanner)?;

    if unique_count == 0 {
        eprintln!("WARNING: No valid minimizers found.");
        return Ok(());
    }

    eprintln!("\nBuilding BF with FPR of {} ...", fpr);

    let bloom = AtomicBloom::new(unique_count, fpr);
    let bloom_ref = &bloom;
    let (seq_tx, seq_rx) = crossbeam_channel::bounded::<Vec<Vec<u8>>>(16);
    let fasta_path = task.fasta_path.clone();

    let reader_handle = std::thread::spawn(move || -> Result<()> {
        let mut reader = needletail::parse_fastx_file(&fasta_path)
            .map_err(|e| anyhow::anyhow!("Cannot open FASTA {}: {}", fasta_path, e))?;

        let mut batch = Vec::with_capacity(1000);
        let mut batch_bytes = 0usize;

        while let Some(record) = reader.next() {
            let rec = record.map_err(|e| anyhow::anyhow!("FASTA error: {}", e))?;
            let seq = rec.seq().into_owned();
            batch_bytes += seq.len();
            batch.push(seq);
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
    eprintln!("Inserted {} total minimizers (including duplicates).", inserted);

    let num_hashes = bloom.num_hashes();
    let num_bits = bloom.num_bits();
    let domain_filter = DomainBloomFilter::new(bloom.into_bytes(), num_hashes, num_bits);
    let out_path = format!("{}.{}", db_prefix, task.output_ext);

    let f = File::create(&out_path)?;
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, f);
    bincode::serialize_into(&mut writer, &domain_filter)?;
    Ok(())
}

fn estimate_unique_minimizers(fasta_path: &str, scanner: &MinimizerScanner) -> Result<(usize, u64)> {
    eprintln!("Estimating cardinality with HLL++...");
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