use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::Mutex;
use std::time::{Duration, Instant};

use anyhow::Result;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::database::{
    delta_vbyte_decode, delta_vbyte_encode, AccessionRegistry, KmerDatabaseBuilder, MinimizerEntry,
};
use crate::minimizer::{create_spaced_seed_mask, MinimizerScanner};
use crate::taxonomy::{load_target_taxids, TargetTaxIDManager, TaxonomyTree};

pub struct BuildConfig {
    pub kraken_file: String,
    pub fasta_file: String,
    pub targets_file: String,
    pub nodes_file: String,
    pub db_prefix: String,
    pub threads: usize,
    pub track_accessions: bool,
    pub k: usize,
    pub l: usize
}

const SPACED_PATTERN: &str = "1111111111111111111110101010101";
const TOGGLE_MASK: u64 = 0xe37e28c4271b5a2d;
const NUM_SHARDS: usize = 64;
const HUMAN_TAXID: u32 = 9606;

const BATCH_SIZE_LIMIT_BYTES: usize = 50 * 1024 * 1024; // 50 MB
const BATCH_SIZE_LIMIT_COUNT: usize = 1000;


pub fn run_build(config: BuildConfig) -> Result<()> {
    let spaced_seed_mask = create_spaced_seed_mask(SPACED_PATTERN);

    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok();

    let total_start = Instant::now();

    eprintln!("Loading Taxonomy...");
    let taxonomy = TaxonomyTree::load(&config.nodes_file)?;

    eprintln!("Loading Target TaxIDs...");
    let targets = load_target_taxids(&config.targets_file)?;
    let taxid_manager = TargetTaxIDManager::new(&targets, &taxonomy);

    let relevant_taxids = taxid_manager.all_relevant_taxids();

    let mut accession_registry = if config.track_accessions {
        eprintln!("\nSorting accessions by taxonomic lineage...");
        let reg = presort_accessions_by_lineage(&config.kraken_file, &taxonomy, &taxid_manager, &relevant_taxids)?;
        Some(reg)
    } else {
        None
    };

    let scanner = MinimizerScanner::new(config.k, config.l, spaced_seed_mask, TOGGLE_MASK);

    let mut minimizer_shards = build_minimizer_shards_streaming(
        &config,
        &taxid_manager,
        &relevant_taxids,
        &mut accession_registry,
        &scanner,
    )?;

    // Filter conflicts directly in shards (parallel)
    eprintln!("\nFiltering conflicts and finalizing database...");
    let total_before: usize = minimizer_shards.iter().map(|m| m.len()).sum();
    
    minimizer_shards.par_iter_mut().for_each(|map| {
        map.retain(|_, entry| !entry.conflicted);
        map.shrink_to_fit();
    });

    let total_after: usize = minimizer_shards.iter().map(|m| m.len()).sum();
    eprintln!(
        "  {} total minimizers, {} conflicted removed, {} remaining",
        total_before,
        total_before - total_after,
        total_after
    );

    eprintln!("Building MPHF and saving...");
    let db = KmerDatabaseBuilder::new(config.k, config.l, spaced_seed_mask, TOGGLE_MASK, config.track_accessions)
        .build_from_minimizers(minimizer_shards,&taxonomy)?;

    db.save(&config.db_prefix)?;

    if let Some(reg) = accession_registry {
        reg.save(&format!("{}.accessions", config.db_prefix))?;
        eprintln!("Saved {} accessions", reg.len());
    }

    eprintln!("Done! Total time: {:.2}s", total_start.elapsed().as_secs_f64());
    Ok(())
}

fn presort_accessions_by_lineage(
    kraken_file: &str,
    taxonomy: &TaxonomyTree,
    taxid_manager: &TargetTaxIDManager,
    relevant_taxids: &HashSet<u32>,
) -> Result<AccessionRegistry> {
    let reader = BufReader::with_capacity(4 * 1024 * 1024, File::open(kraken_file)?);

    // Collect (seq_id, kraken_taxid) for sequences needing accession tracking.
    // We only need to know: does this sequence have any non-human target hit?
    // And what is its kraken2 classification taxid (column 3)?
    let mut accession_seqs: Vec<(String, u32)> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut line_count = 0u64;

    for line_res in reader.lines() {
        let line = line_res?;
        line_count += 1;

        if line_count % 1_000_000 == 0 {
            eprint!("\r  Pre-scan: {}M lines...", line_count / 1_000_000);
        }

        let mut parts = line.split('\t');
        parts.next(); // skip classification status
        let seq_id = match parts.next() { Some(v) => v, None => continue };
        let kraken_taxid: u32 = match parts.next().and_then(|s| s.trim().parse().ok()) {
            Some(v) => v,
            None => continue,
        };
        if parts.next().is_none() { continue; } 
        let lca_mapping = parts.next().unwrap_or("");

        // Check if this sequence has any non-human target hits
        let mut has_non_human_target = false;
        for part in lca_mapping.split_whitespace() {
            let colon_pos = match part.find(':') { Some(p) => p, None => continue };
            let taxid_str = &part[..colon_pos];
            if taxid_str == "cov" { continue; }

            let taxid: u32 = match taxid_str.parse() { Ok(v) => v, Err(_) => continue };

            if relevant_taxids.contains(&taxid) {
                if let Some(target_taxid) = taxid_manager.get_target(taxid) {
                    if target_taxid != HUMAN_TAXID {
                        has_non_human_target = true;
                        break;
                    }
                }
            }
        }

        if has_non_human_target && seen.insert(seq_id.to_string()) {
            accession_seqs.push((seq_id.to_string(), kraken_taxid));
        }
    }

    eprintln!("\r  Pre-scan: {} lines, {} accessions to sort", line_count, accession_seqs.len());

    // Cache lineage paths per unique taxid
    let mut lineage_cache: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for &(_, taxid) in &accession_seqs {
        lineage_cache
            .entry(taxid)
            .or_insert_with(|| taxonomy.lineage_path(taxid));
    }

    // Stable sort by lineage path — preserves encounter order within same taxid
    accession_seqs.sort_by(|a, b| {
        lineage_cache[&a.1].cmp(&lineage_cache[&b.1])
    });

    // Pre-assign IDs in sorted order
    let mut registry = AccessionRegistry::new();
    for (seq_id, _) in &accession_seqs {
        registry.get_or_create(seq_id);
    }

    eprintln!("  Pre-assigned {} accession IDs in taxonomic order", registry.len());
    Ok(registry)
}

fn build_minimizer_shards_streaming(
    config: &BuildConfig,
    taxid_manager: &TargetTaxIDManager,
    relevant_taxids: &HashSet<u32>,
    accession_registry: &mut Option<AccessionRegistry>,
    scanner: &MinimizerScanner,
) -> Result<Vec<FxHashMap<u64, MinimizerEntry>>> {
    
    eprintln!("\nPre-scanning Kraken output for targets...");
    
    // Key: Sequence ID -> Value: (Optional Accession ID, Vec<(start_idx, count, target_taxid)>)
    let mut job_map: FxHashMap<String, (Option<u32>, Vec<(usize, usize, u32)>)> = FxHashMap::default();
    let kraken_reader = BufReader::with_capacity(4 * 1024 * 1024, File::open(&config.kraken_file)?);
    
    let mut total_kraken_lines = 0;
    let mut relevant_seq_count = 0;

    for line_res in kraken_reader.lines() {
        let line = line_res?;
        total_kraken_lines += 1;
        
        if total_kraken_lines % 1_000_000 == 0 {
            eprint!("\r  Scanned {}M lines...", total_kraken_lines / 1_000_000);
        }

        let mut parts = line.split('\t');
        parts.next(); 
        let seq_id = match parts.next() { Some(v) => v, None => continue };
        if parts.next().is_none() || parts.next().is_none() { continue; }
        let lca_mapping = parts.next().unwrap_or("");

        let mut kmer_read_index: usize = 0;
        let mut jobs = Vec::new();

        for part in lca_mapping.split_whitespace() {
            let colon_pos = match part.find(':') { Some(p) => p, None => continue };
            let taxid_str = &part[..colon_pos];
            let count_str = &part[colon_pos + 1..];

            if taxid_str == "cov" { continue; }

            let count: usize = match count_str.parse() { Ok(v) => v, Err(_) => continue };
            let taxid: u32 = match taxid_str.parse() {
                Ok(v) => v,
                Err(_) => { kmer_read_index += count; continue; }
            };

            if relevant_taxids.contains(&taxid) {
                if let Some(target_taxid) = taxid_manager.get_target(taxid) {
                    jobs.push((kmer_read_index, count, target_taxid));
                }
            }
            kmer_read_index += count;
        }

        if !jobs.is_empty() {
             let needs_accession = jobs.iter().any(|&(_, _, tid)| tid != HUMAN_TAXID);

             if let std::collections::hash_map::Entry::Vacant(e) = job_map.entry(seq_id.to_string()) {
                // If pre-sorted registry exists, get_or_create just returns the
                // already-assigned ID. Otherwise assigns in encounter order.
                let acc_id = if needs_accession {
                    if let Some(ref mut reg) = accession_registry {
                        Some(reg.get_or_create(seq_id))
                    } else { None }
                } else { None };
                e.insert((acc_id, jobs));
                relevant_seq_count += 1;
            } else {
                let entry = job_map.get_mut(seq_id).unwrap();
                entry.1.extend(jobs);
                if entry.0.is_none() && needs_accession {
                    if let Some(ref mut reg) = accession_registry {
                        entry.0 = Some(reg.get_or_create(seq_id));
                    }
                }
            }
        }
    }
    eprintln!("\r  Scanned {} lines. Found {} relevant sequences.", total_kraken_lines, relevant_seq_count);
    eprintln!("\nStreaming FASTA and extracting minimizers...");

    // Job payload sent through channel — producer extracts from job_map, no mutex needed by consumers
    type JobPayload = Vec<(Vec<u8>, Option<u32>, Vec<(usize, usize, u32)>)>;
    
    // Bounded channel
    let (tx, rx) = crossbeam_channel::bounded::<JobPayload>(32);
    let fasta_file = config.fasta_file.clone();
    
    // Producer Thread — owns the job_map, no Arc/Mutex needed
    let reader_handle = std::thread::spawn(move || -> Result<()> {
        let mut reader = needletail::parse_fastx_file(&fasta_file)
            .map_err(|e| anyhow::anyhow!("Cannot open FASTA file: {}", e))?;
        
        let mut batch: JobPayload = Vec::new();
        let mut batch_bytes = 0;
        let mut total_bytes_processed = 0u64;
        let mut last_log = Instant::now();

        while let Some(record) = reader.next() {
            let rec = record.map_err(|e| anyhow::anyhow!("FASTA record error: {}", e))?;
            let id_full = rec.id();
            let id_str = std::str::from_utf8(id_full).unwrap_or("");
            let seq_id = id_str.split_whitespace().next().unwrap_or("");

            // Extract AND REMOVE the job here in the producer — consumers never touch job_map
            if let Some((acc_id, jobs)) = job_map.remove(seq_id) {
                let seq_data = rec.seq().into_owned();
                batch_bytes += seq_data.len();
                batch.push((seq_data, acc_id, jobs));

                if batch.len() >= BATCH_SIZE_LIMIT_COUNT || batch_bytes >= BATCH_SIZE_LIMIT_BYTES {
                    if tx.send(batch).is_err() { return Ok(()); }
                    batch = Vec::new();
                    batch_bytes = 0;
                }
            }

            total_bytes_processed += rec.num_bases() as u64; 
            if last_log.elapsed() > Duration::from_secs(5) {
                eprint!("\r  Processed {:.2} GB...", total_bytes_processed as f64 / 1_073_741_824.0);
                last_log = Instant::now();
            }

            // Early exit: all jobs consumed
            if job_map.is_empty() {
                break;
            }
        }
        
        if !batch.is_empty() {
            let _ = tx.send(batch);
        }
        eprintln!("\r  Finished reading FASTA ({:.2} GB). Waiting for workers...", total_bytes_processed as f64 / 1_073_741_824.0);
        Ok(())
    });

    let global = ShardedMap::new();
    let track_accessions = config.track_accessions;
    let k = scanner.k();

    // Consumer Threads (Rayon) — no mutex, just process the payloads
    rx.into_iter().par_bridge().for_each(|batch| {
        let mut local_map: FxHashMap<u64, RawLocalEntry> = FxHashMap::default();
        let mut minimizer_buf: Vec<u64> = Vec::new();

        for (seq_bytes, acc_id_opt, jobs) in &batch {
            let seq_len = seq_bytes.len();

            for &(kmer_read_index, count, target_taxid) in jobs {
                let region_start = kmer_read_index;
                let region_end = (kmer_read_index + count - 1 + k).min(seq_len);

                if region_start < seq_len && region_end <= seq_len && region_end > region_start {
                    scanner.scan_into(&seq_bytes[region_start..region_end], &mut minimizer_buf);

                    if !minimizer_buf.is_empty() {
                        let aid = if track_accessions && target_taxid != HUMAN_TAXID {
                            *acc_id_opt
                        } else {
                            None
                        };

                        for &m in &minimizer_buf {
                            let entry = local_map.entry(m).or_insert_with(|| RawLocalEntry {
                                taxid: target_taxid,
                                conflicted: false,
                                accession_ids: Vec::new(),
                            });
                            if entry.taxid != target_taxid {
                                entry.conflicted = true;
                            }
                            if let Some(id) = aid {
                                entry.accession_ids.push(id);
                            }
                        }
                    }
                }
            }
        }

        if track_accessions {
            for entry in local_map.values_mut() {
                entry.accession_ids.sort_unstable();
                entry.accession_ids.dedup();
            }
        }

        let mut by_shard: Vec<Vec<(u64, RawLocalEntry)>> = (0..NUM_SHARDS).map(|_| Vec::new()).collect();
        for (minimizer, entry) in local_map {
            by_shard[shard_for(minimizer)].push((minimizer, entry));
        }

        let mut merge_buf: Vec<u32> = Vec::new();
        for (shard_idx, entries) in by_shard.into_iter().enumerate() {
            if entries.is_empty() { continue; }
            let mut shard = global.shards[shard_idx].lock().unwrap();
            for (minimizer, local) in entries {
                merge_into_shard(&mut shard, minimizer, local, &mut merge_buf, track_accessions);
            }
        }
    });

    match reader_handle.join() {
        Ok(Ok(())) => {},
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("FASTA reader thread panicked"),
    }

    Ok(global.into_shards())
}

#[inline(always)]
fn shard_for(minimizer: u64) -> usize {
    (minimizer >> 58) as usize & (NUM_SHARDS - 1)
}

struct ShardedMap {
    shards: Vec<Mutex<FxHashMap<u64, MinimizerEntry>>>,
}

impl ShardedMap {
    fn new() -> Self {
        Self {
            shards: (0..NUM_SHARDS).map(|_| Mutex::new(FxHashMap::default())).collect(),
        }
    }

    fn into_shards(self) -> Vec<FxHashMap<u64, MinimizerEntry>> {
        self.shards.into_iter().map(|m| m.into_inner().unwrap()).collect()
    }
}

struct RawLocalEntry {
    taxid: u32,
    conflicted: bool,
    accession_ids: Vec<u32>,
}

fn merge_into_shard(
    shard: &mut FxHashMap<u64, MinimizerEntry>,
    minimizer: u64,
    local: RawLocalEntry,
    merge_buf: &mut Vec<u32>,
    track_accessions: bool,
) {
    let global = shard.entry(minimizer).or_insert_with(|| MinimizerEntry {
        taxid: local.taxid,
        conflicted: local.conflicted,
        accessions_vbyte: Box::from([]),
    });

    if global.taxid != local.taxid || local.conflicted {
        global.conflicted = true;
    }

    if !track_accessions || local.accession_ids.is_empty() {
        return;
    }

    if global.accessions_vbyte.is_empty() {
        global.accessions_vbyte = delta_vbyte_encode(&local.accession_ids).into_boxed_slice();
        return;
    }

    let mut pos = 0;
    let a = delta_vbyte_decode(&global.accessions_vbyte, &mut pos);
    let b = &local.accession_ids;

    merge_buf.clear();
    merge_buf.reserve(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        use std::cmp::Ordering;
        match a[i].cmp(&b[j]) {
            Ordering::Less    => { merge_buf.push(a[i]); i += 1; }
            Ordering::Greater => { merge_buf.push(b[j]); j += 1; }
            Ordering::Equal   => { merge_buf.push(a[i]); i += 1; j += 1; }
        }
    }
    merge_buf.extend_from_slice(&a[i..]);
    merge_buf.extend_from_slice(&b[j..]);

    global.accessions_vbyte = delta_vbyte_encode(merge_buf).into_boxed_slice();
}