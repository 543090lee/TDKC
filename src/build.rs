use std::collections::{BTreeSet, HashMap, HashSet};
use std::sync::{Arc, RwLock};
use std::time::{Duration, Instant};

use anyhow::Result;
use boomphf::Mphf;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::database::{AccessionRegistry, EqClassAccessions, KmerDatabase};
use crate::minimizer::{create_spaced_seed_mask, MinimizerScanner};
use crate::prep::load_prelim_map;
use crate::taxonomy::{load_target_taxids, BfsTaxonomy, TargetTaxIDManager, TaxonomyTree};

pub struct BuildConfig {
    pub fasta_file: String,
    pub target_fasta_file: String,
    pub prelim_map_file: String,
    pub targets_file: String,
    pub nodes_file: String,
    pub db_prefix: String,
    pub threads: usize,
    pub track_accessions: bool,
    pub k: usize,
    pub l: usize,
}

const SPACED_PATTERN: &str = "1111111111111111101010101010101";
const TOGGLE_MASK: u64 = 0xe37e28c4271b5a2d;

const NUM_SHARDS: usize = 256;
const SHARD_SHIFT: usize = 56;

const BATCH_ITEM_LIMIT: usize = 2048;
const BATCH_BYTE_LIMIT: usize = 256 * 1024 * 1024;

const PAR_CHUNK_SIZE: usize = 64;

const SEG_TARGET_LEN: usize = 4 * 1024 * 1024;

pub fn run_build(config: BuildConfig) -> Result<()> {
    let spaced_seed_mask = create_spaced_seed_mask(SPACED_PATTERN);

    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok();

    let total_start = Instant::now();

    eprintln!("Phase 0: Loading taxonomy...");
    let taxonomy = TaxonomyTree::load(&config.nodes_file)?;
    log_mem("Loaded TaxonomyTree");

    eprintln!("Loading target TaxIDs...");
    let targets = load_target_taxids(&config.targets_file)?;
    let taxid_manager = TargetTaxIDManager::new(&targets, &taxonomy);
    let relevant_taxids = taxid_manager.all_relevant_taxids();
    log_mem("Loaded Targets and Manager");

    eprintln!("Building BFS taxonomy with internal IDs...");
    let bfs_tax = BfsTaxonomy::build(&taxonomy, &relevant_taxids);
    log_mem("Built BFS Taxonomy");

    eprintln!("\nLoading prelim_map...");
    let acc_to_taxid = load_prelim_map(&config.prelim_map_file)?;
    log_mem("Loaded prelim_map");

    let acc_to_internal: HashMap<String, u32> = acc_to_taxid
        .iter()
        .filter_map(|(acc, ext_taxid): (&String, &u32)| {
            bfs_tax.to_internal(*ext_taxid).map(|int_id| (acc.clone(), int_id))
        })
        .collect();
    eprintln!(
        "  {} accessions have taxids in taxonomy tree",
        acc_to_internal.len()
    );
    log_mem("Built acc_to_internal mapping");

    let target_accessions: HashSet<String> = acc_to_internal
        .iter()
        .filter(|(_, int_id): &(&String, &u32)| bfs_tax.is_relevant(**int_id))
        .map(|(acc, _): (&String, &u32)| acc.clone())
        .collect();
    eprintln!(
        "  {} accessions belong to target clades",
        target_accessions.len()
    );
    log_mem("Filtered target_accessions");

    let scanner = MinimizerScanner::new(config.k, config.l, spaced_seed_mask, TOGGLE_MASK);

    let mut accession_registry = if config.track_accessions {
        eprintln!("\nSorting accessions by taxonomic lineage...");
        let reg = presort_accessions_by_lineage(
            &target_accessions,
            &acc_to_taxid,
            &taxonomy,
            &taxid_manager,
        )?;
        log_mem("Built AccessionRegistry");
        Some(reg)
    } else {
        None
    };

    eprintln!("\nPhase 2: Extracting minimizers from target sub-FASTA...");
    let global = ShardedMinimizerMap::new();
    log_mem("Initialized ShardedMinimizerMap (Empty)");

    phase2_extract_target_minimizers(
        &config.target_fasta_file,
        &target_accessions,
        &acc_to_internal,
        &bfs_tax,
        &scanner,
        &global,
        config.track_accessions,
        &mut accession_registry,
    )?;
    log_mem("Finished Phase 2 (Peak minimizer extraction)");

    let total_target_minimizers: usize = global
        .shards
        .iter()
        .map(|s| s.read().unwrap().len())
        .sum();
    eprintln!("  {} base target minimizers", total_target_minimizers);

    eprintln!("\nPhase 3: Challenging minimizers with background sequences...");
    let to_remove = phase3_collect_background_hits(
        &config.fasta_file,
        &target_accessions,
        &scanner,
        &global,
    )?;

    let remove_count = to_remove.len();
    eprintln!("  {} minimizers found in background, removing...", remove_count);
    log_mem("Collected Phase 3 background hits");

    phase3_bulk_remove(&global, &to_remove);
    drop(to_remove);
    log_mem("Completed Phase 3 bulk remove");

    let total_surviving: usize = global
        .shards
        .iter()
        .map(|s| s.read().unwrap().len())
        .sum();
    eprintln!(
        "  {} minimizers survived ({} removed by challenge)",
        total_surviving,
        total_target_minimizers.saturating_sub(total_surviving)
    );

    eprintln!("\nPhase 4: Building database directly from sharded map...");
    let db = build_database_from_shards(
        global,
        &bfs_tax,
        &taxid_manager,
        &taxonomy,
        config.k,
        config.l,
        spaced_seed_mask,
        TOGGLE_MASK,
        config.track_accessions,
    )?;
    log_mem("Built Final KmerDatabase in Memory");

    db.save(&config.db_prefix)?;
    log_mem("Saved KmerDatabase to disk");

    if let Some(reg) = accession_registry {
        reg.save(&format!("{}.accessions", config.db_prefix))?;
        eprintln!("Saved {} accessions", reg.len());
    }

    eprintln!(
        "\nDone! Total time: {:.2}s",
        total_start.elapsed().as_secs_f64()
    );
    Ok(())
}

fn build_database_from_shards(
    global: ShardedMinimizerMap,
    bfs_tax: &BfsTaxonomy,
    taxid_manager: &TargetTaxIDManager,
    taxonomy: &TaxonomyTree,
    k: usize,
    l: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
    track_accessions: bool,
) -> Result<KmerDatabase> {
    let start = Instant::now();
    eprintln!("\nBuilding k-mer database");

    // ── Pass 1: Scan shards to collect unique taxids and count ───────────────
    eprintln!("\nCreating TaxID mapping (scanning shards)...");
    let mut unique_taxids = BTreeSet::new();
    let mut num_minimizers: usize = 0;

    for shard_lock in &global.shards {
        let shard = shard_lock.read().unwrap();
        for entry in shard.values() {
            if !bfs_tax.is_relevant(entry.int_taxid) {
                continue;
            }
            let ext_taxid = bfs_tax.to_external(entry.int_taxid);
            let report_taxid = taxid_manager.get_target(ext_taxid).unwrap_or(ext_taxid);
            unique_taxids.insert(report_taxid);
            num_minimizers += 1;
        }
    }

    eprintln!("  Found {} unique TaxIDs", unique_taxids.len());
    eprintln!("  {} unique minimizers", num_minimizers);

    if unique_taxids.len() > 255 {
        anyhow::bail!(
            "Too many unique taxids ({}) for u8 storage (max 255)",
            unique_taxids.len()
        );
    }
    if num_minimizers == 0 {
        anyhow::bail!("No minimizers found");
    }

    let mut index_to_taxid: Vec<u32> = Vec::new();
    let mut taxid_to_index: HashMap<u32, u8> = HashMap::new();
    for (i, &taxid) in unique_taxids.iter().enumerate() {
        index_to_taxid.push(taxid);
        taxid_to_index.insert(taxid, i as u8);
    }

    eprintln!("\nGenerating Ancestry Matrix...");
    let n = index_to_taxid.len();
    let mut ancestor_matrix = vec![0u8; n * n];
    for i in 0..n {
        for j in 0..n {
            let parent_taxid = index_to_taxid[i];
            let child_taxid = index_to_taxid[j];
            if parent_taxid == child_taxid {
                ancestor_matrix[i * n + j] = 1;
            } else {
                let path = taxonomy.lineage_path(child_taxid);
                if path.contains(&parent_taxid) {
                    ancestor_matrix[i * n + j] = 1;
                }
            }
        }
    }

    eprintln!("\nCollecting keys for MPHF ({} keys)...", num_minimizers);
    let mut keys: Vec<u64> = Vec::with_capacity(num_minimizers);
    for shard_lock in &global.shards {
        let shard = shard_lock.read().unwrap();
        for (&minimizer, entry) in shard.iter() {
            if bfs_tax.is_relevant(entry.int_taxid) {
                keys.push(minimizer);
            }
        }
    }
    log_mem("Collected MPHF keys");

    eprintln!("Building MPHF...");
    let mphf = Mphf::new(2.0, &keys);
    drop(keys);
    log_mem("Built MPHF (keys dropped)");

    eprintln!("\nPopulating arrays (consuming shards to free memory)...");
    let mut fingerprints = vec![0u16; num_minimizers];
    let mut taxid_indices = vec![0u8; num_minimizers];
    let mut accessions = track_accessions
        .then(|| EqClassAccessions::new_empty(num_minimizers));
    let mut hash_to_class: FxHashMap<u64, u32> = FxHashMap::default();
    let mut build_map: FxHashMap<u32, u32> = FxHashMap::default();

    for (shard_idx, shard_lock) in global.shards.into_iter().enumerate() {
        let shard = shard_lock.into_inner().unwrap();

        for (minimizer, entry) in shard {
            if !bfs_tax.is_relevant(entry.int_taxid) {
                continue;
            }

            let ext_taxid = bfs_tax.to_external(entry.int_taxid);
            let report_taxid = taxid_manager.get_target(ext_taxid).unwrap_or(ext_taxid);

            let idx = mphf.hash(&minimizer) as usize;
            if idx >= num_minimizers {
                continue;
            }

            fingerprints[idx] = compute_fingerprint(minimizer);
            taxid_indices[idx] = taxid_to_index[&report_taxid];

            if let Some(ref mut acc) = accessions {
                if let Some(mut boxed_vec) = entry.accession_ids {
                    if !boxed_vec.is_empty() {
                        boxed_vec.sort_unstable();
                        boxed_vec.dedup();
                        acc.add_accessions(idx, &boxed_vec, &mut hash_to_class, &mut build_map);
                    }
                }
            }
        }

        if (shard_idx + 1) % 64 == 0 {
            log_mem(&format!("  Consumed {}/{} shards", shard_idx + 1, NUM_SHARDS));
        }
    }
    log_mem("All shards consumed, arrays populated");

    drop(hash_to_class);

    if let Some(ref mut acc) = accessions {
        eprintln!(
            "  Equivalence classes: {} unique sets for accessions",
            acc.num_classes - 1
        );
        acc.finalize_for_save(build_map);
    } else {
        drop(build_map);
    }

    eprintln!("\nDone building in {:.2}s", start.elapsed().as_secs_f64());

    Ok(KmerDatabase::new(
        k,
        l,
        spaced_seed_mask,
        toggle_mask,
        num_minimizers,
        mphf,
        fingerprints,
        taxid_indices,
        index_to_taxid,
        ancestor_matrix,
        accessions,
    ))
}

#[inline]
fn compute_fingerprint(mut kmer: u64) -> u16 {
    kmer ^= kmer >> 33;
    kmer = kmer.wrapping_mul(0xff51afd7ed558ccd);
    kmer ^= kmer >> 33;
    (kmer & 0xFFFF) as u16
}


fn presort_accessions_by_lineage(
    target_accessions: &HashSet<String>,
    acc_to_taxid: &HashMap<String, u32>,
    taxonomy: &TaxonomyTree,
    taxid_manager: &TargetTaxIDManager,
) -> Result<AccessionRegistry> {
    let mut accession_seqs: Vec<(String, u32)> = Vec::new();

    for acc in target_accessions {
        if let Some(&ext_taxid) = acc_to_taxid.get(acc) {
            let report_taxid = taxid_manager.get_target(ext_taxid).unwrap_or(ext_taxid);
            if report_taxid != 9606 {
                accession_seqs.push((acc.clone(), report_taxid));
            }
        }
    }

    let mut lineage_cache: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for &(_, taxid) in &accession_seqs {
        lineage_cache
            .entry(taxid)
            .or_insert_with(|| taxonomy.lineage_path(taxid));
    }

    accession_seqs.sort_by(|a, b| lineage_cache[&a.1].cmp(&lineage_cache[&b.1]));

    let mut registry = AccessionRegistry::new();
    for (seq_id, _) in &accession_seqs {
        registry.get_or_create(seq_id);
    }

    eprintln!(
        "  Pre-assigned {} accession IDs in taxonomic order",
        registry.len()
    );
    Ok(registry)
}


struct BuildEntry {
    int_taxid: u32,
    accession_ids: Option<Box<Vec<u32>>>,
}

struct ShardedMinimizerMap {
    shards: Vec<RwLock<FxHashMap<u64, BuildEntry>>>,
}

impl ShardedMinimizerMap {
    fn new() -> Self {
        Self {
            shards: (0..NUM_SHARDS)
                .map(|_| RwLock::new(FxHashMap::default()))
                .collect(),
        }
    }
}

#[inline(always)]
fn shard_for(minimizer: u64) -> usize {
    (minimizer >> SHARD_SHIFT) as usize & (NUM_SHARDS - 1)
}

struct LocalShardAccum {
    shards: Vec<Vec<(u64, u32)>>,
}

impl LocalShardAccum {
    fn new() -> Self {
        Self {
            shards: (0..NUM_SHARDS).map(|_| Vec::new()).collect(),
        }
    }

    #[inline]
    fn push(&mut self, minimizer: u64, int_taxid: u32) {
        self.shards[shard_for(minimizer)].push((minimizer, int_taxid));
    }

    fn flush_insert(
        &mut self,
        global: &ShardedMinimizerMap,
        bfs_tax: &BfsTaxonomy,
    ) {
        for (shard_idx, entries) in self.shards.iter_mut().enumerate() {
            if entries.is_empty() {
                continue;
            }

            entries.sort_unstable_by_key(|&(m, _)| m);
            dedup_lca(entries, bfs_tax);

            let mut shard = global.shards[shard_idx].write().unwrap();
            for &(minimizer, int_taxid) in entries.iter() {
                match shard.get_mut(&minimizer) {
                    Some(existing) => {
                        let lca = bfs_tax.lca(existing.int_taxid, int_taxid);
                        existing.int_taxid = lca;
                        if !bfs_tax.is_relevant(lca) {
                            existing.accession_ids = None;
                        }
                    }
                    None => {
                        shard.insert(
                            minimizer,
                            BuildEntry {
                                int_taxid,
                                accession_ids: None,
                            },
                        );
                    }
                }
            }
            entries.clear();
        }
    }
}

#[inline]
fn dedup_lca(entries: &mut Vec<(u64, u32)>, bfs_tax: &BfsTaxonomy) {
    entries.dedup_by(|a, b| {
        if a.0 == b.0 {
            b.1 = bfs_tax.lca(a.1, b.1);
            true
        } else {
            false
        }
    });
}

// ─── Work items ─────────────────────────────────────────────────────────────

struct WorkItem {
    seq_data: Arc<Vec<u8>>,
    start: usize,
    end: usize,
    int_taxid: u32,
}

impl WorkItem {
    #[inline]
    fn seq_slice(&self) -> &[u8] {
        &self.seq_data[self.start..self.end]
    }
}

struct AccRecord {
    seq_data: Arc<Vec<u8>>,
    acc_name: String,
}

fn segment_ranges(seq_len: usize, k: usize) -> Vec<(usize, usize)> {
    if seq_len <= SEG_TARGET_LEN {
        return vec![(0, seq_len)];
    }

    let overlap = k.saturating_sub(1);
    let stride = SEG_TARGET_LEN - overlap;
    let mut ranges = Vec::new();
    let mut start = 0;

    while start < seq_len {
        let end = (start + SEG_TARGET_LEN).min(seq_len);
        ranges.push((start, end));
        if end == seq_len {
            break;
        }
        start += stride;
    }

    ranges
}


fn phase2_extract_target_minimizers(
    target_fasta_path: &str,
    target_accessions: &HashSet<String>,
    acc_to_internal: &HashMap<String, u32>,
    bfs_tax: &BfsTaxonomy,
    scanner: &MinimizerScanner,
    global: &ShardedMinimizerMap,
    track_accessions: bool,
    accession_registry: &mut Option<AccessionRegistry>,
) -> Result<()> {
    let k = scanner.k();

    let (tx, rx) = crossbeam_channel::bounded::<(Vec<WorkItem>, Vec<AccRecord>)>(4);
    let fasta_file = target_fasta_path.to_string();
    let acc_int = acc_to_internal.clone();
    let target_acc = target_accessions.clone();
    let want_accessions = track_accessions;
    let human_int_taxid = bfs_tax.to_internal(9606);

    let reader_handle = std::thread::spawn(move || -> Result<()> {
        let mut reader = needletail::parse_fastx_file(&fasta_file)
            .map_err(|e| anyhow::anyhow!("Cannot open target FASTA: {}", e))?;

        let mut work_items: Vec<WorkItem> = Vec::new();
        let mut acc_records: Vec<AccRecord> = Vec::new();
        let mut batch_bytes = 0usize;
        let mut total_seqs = 0u64;
        let mut total_segments = 0u64;
        let mut last_log = Instant::now();

        while let Some(record) = reader.next() {
            let rec = record.map_err(|e| anyhow::anyhow!("FASTA error: {}", e))?;
            let id_full = std::str::from_utf8(rec.id()).unwrap_or("");
            let accession = id_full.split_whitespace().next().unwrap_or("");
            total_seqs += 1;

            let int_taxid = match acc_int.get(accession) {
                Some(&id) => id,
                None => continue,
            };

            let seq_data = Arc::new(rec.seq().into_owned());
            let seq_len = seq_data.len();

            if want_accessions && target_acc.contains(accession) {
                if Some(int_taxid) != human_int_taxid {
                    acc_records.push(AccRecord {
                        seq_data: Arc::clone(&seq_data),
                        acc_name: accession.to_string(),
                    });
                }
            }

            let ranges = segment_ranges(seq_len, k);
            let num_segs = ranges.len();
            total_segments += num_segs as u64;

            for (start, end) in ranges {
                batch_bytes += end - start;
                work_items.push(WorkItem {
                    seq_data: Arc::clone(&seq_data),
                    start,
                    end,
                    int_taxid,
                });
            }

            if num_segs > 1 && seq_len > 100_000_000 {
                eprintln!(
                    "    Split {} ({:.0} MB) into {} segments",
                    accession,
                    seq_len as f64 / 1_048_576.0,
                    num_segs
                );
            }

            if work_items.len() >= BATCH_ITEM_LIMIT || batch_bytes >= BATCH_BYTE_LIMIT {
                if tx
                    .send((
                        std::mem::replace(&mut work_items, Vec::new()),
                        std::mem::replace(&mut acc_records, Vec::new()),
                    ))
                    .is_err()
                {
                    return Ok(());
                }
                batch_bytes = 0;
            }

            if last_log.elapsed() > Duration::from_secs(5) {
                eprint!(
                    "\r  Phase 2: {} seqs ({} work items)...",
                    total_seqs, total_segments
                );
                last_log = Instant::now();
            }
        }

        if !work_items.is_empty() || !acc_records.is_empty() {
            let _ = tx.send((work_items, acc_records));
        }
        eprintln!(
            "\r  Phase 2: {} target seqs, {} work items",
            total_seqs, total_segments
        );
        Ok(())
    });

    for (work_items, acc_records) in rx {
        work_items.par_chunks(PAR_CHUNK_SIZE).for_each(|chunk| {
            let mut accum = LocalShardAccum::new();
            let mut minimizer_buf: Vec<u64> = Vec::new();

            for item in chunk {
                scanner.scan_into(item.seq_slice(), &mut minimizer_buf);
                for &m in &minimizer_buf {
                    if m != u64::MAX {
                        accum.push(m, item.int_taxid);
                    }
                }
            }

            accum.flush_insert(global, bfs_tax);
        });

        if track_accessions {
            if let Some(ref mut registry) = accession_registry {
                let mut minimizer_buf: Vec<u64> = Vec::new();

                for arec in &acc_records {
                    let acc_id = registry.get_or_create(&arec.acc_name);
                    scanner.scan_into(&arec.seq_data, &mut minimizer_buf);

                    let mut by_shard: Vec<Vec<u64>> =
                        (0..NUM_SHARDS).map(|_| Vec::new()).collect();
                    for &m in &minimizer_buf {
                        if m != u64::MAX {
                            by_shard[shard_for(m)].push(m);
                        }
                    }

                    for (shard_idx, mut mins) in by_shard.into_iter().enumerate() {
                        if mins.is_empty() {
                            continue;
                        }

                        mins.sort_unstable();
                        mins.dedup();

                        let mut shard = global.shards[shard_idx].write().unwrap();
                        for m in mins {
                            if let Some(entry) = shard.get_mut(&m) {
                                if bfs_tax.is_relevant(entry.int_taxid) {
                                    entry
                                        .accession_ids
                                        .get_or_insert_with(|| Box::new(Vec::new()))
                                        .push(acc_id);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    match reader_handle.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("Phase 2 reader thread panicked"),
    }

    Ok(())
}
struct Phase3Segment {
    seq_data: Arc<Vec<u8>>,
    start: usize,
    end: usize,
}

impl Phase3Segment {
    #[inline]
    fn seq_slice(&self) -> &[u8] {
        &self.seq_data[self.start..self.end]
    }
}

fn phase3_collect_background_hits(
    fasta_path: &str,
    target_accessions: &HashSet<String>,
    scanner: &MinimizerScanner,
    global: &ShardedMinimizerMap,
) -> Result<Vec<u64>> {
    let k = scanner.k();

    let (tx, rx) = crossbeam_channel::bounded::<Vec<Phase3Segment>>(4);
    let fasta_file = fasta_path.to_string();
    let target_acc = target_accessions.clone();

    let reader_handle = std::thread::spawn(move || -> Result<()> {
        let mut reader = needletail::parse_fastx_file(&fasta_file)
            .map_err(|e| anyhow::anyhow!("Cannot open FASTA: {}", e))?;

        let mut work_items: Vec<Phase3Segment> = Vec::new();
        let mut batch_bytes = 0usize;
        let mut total_seqs = 0u64;
        let mut skipped_target = 0u64;
        let mut background = 0u64;
        let mut last_log = Instant::now();

        while let Some(record) = reader.next() {
            let rec = record.map_err(|e| anyhow::anyhow!("FASTA error: {}", e))?;
            let id_full = std::str::from_utf8(rec.id()).unwrap_or("");
            let accession = id_full.split_whitespace().next().unwrap_or("");
            total_seqs += 1;

            if target_acc.contains(accession) {
                skipped_target += 1;
                if last_log.elapsed() > Duration::from_secs(5) {
                    eprint!(
                        "\r  Phase 3: {}M seqs, {}M background, {}k target skipped...",
                        total_seqs,
                        background,
                        skipped_target
                    );
                    last_log = Instant::now();
                }
                continue;
            }

            background += 1;
            let seq_data = Arc::new(rec.seq().into_owned());
            let seq_len = seq_data.len();

            let ranges = segment_ranges(seq_len, k);
            for (start, end) in ranges {
                batch_bytes += end - start;
                work_items.push(Phase3Segment {
                    seq_data: Arc::clone(&seq_data),
                    start,
                    end,
                });
            }

            if work_items.len() >= BATCH_ITEM_LIMIT || batch_bytes >= BATCH_BYTE_LIMIT {
                if tx
                    .send(std::mem::replace(&mut work_items, Vec::new()))
                    .is_err()
                {
                    return Ok(());
                }
                batch_bytes = 0;
            }

            if last_log.elapsed() > Duration::from_secs(5) {
                eprint!(
                    "\r  Phase 3: {}M seqs, {}M background, {}k target skipped...",
                    total_seqs,
                    background,
                    skipped_target
                );
                last_log = Instant::now();
            }
        }

        if !work_items.is_empty() {
            let _ = tx.send(work_items);
        }
        eprintln!(
            "\r  Phase 3: {} total seqs, {} background scanned, {} target skipped",
            total_seqs, background, skipped_target
        );
        Ok(())
    });

    let mut all_hits: Vec<u64> = Vec::new();

    for work_items in rx {
        let batch_hits: Vec<u64> = work_items
            .par_chunks(PAR_CHUNK_SIZE)
            .fold(
                || Vec::new(),
                |mut local_hits: Vec<u64>, chunk: &[Phase3Segment]| {
                    let mut minimizer_buf: Vec<u64> = Vec::new();

                    for item in chunk {
                        scanner.scan_into(item.seq_slice(), &mut minimizer_buf);
                        for &m in &minimizer_buf {
                            if m != u64::MAX {
                                let shard = global.shards[shard_for(m)].read().unwrap();
                                if shard.contains_key(&m) {
                                    local_hits.push(m);
                                }
                            }
                        }
                    }

                    local_hits
                },
            )
            .reduce(
                || Vec::new(),
                |mut a, mut b| {
                    if a.len() >= b.len() {
                        a.append(&mut b);
                        a
                    } else {
                        b.append(&mut a);
                        b
                    }
                },
            );

        all_hits.extend(batch_hits);
    }

    match reader_handle.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("reader thread panicked"),
    }

    all_hits.sort_unstable();
    all_hits.dedup();

    Ok(all_hits)
}

fn phase3_bulk_remove(global: &ShardedMinimizerMap, to_remove: &[u64]) {
    let mut by_shard: Vec<Vec<u64>> = (0..NUM_SHARDS).map(|_| Vec::new()).collect();
    for &m in to_remove {
        by_shard[shard_for(m)].push(m);
    }

    for (shard_idx, removals) in by_shard.into_iter().enumerate() {
        if removals.is_empty() {
            continue;
        }
        let mut shard = global.shards[shard_idx].write().unwrap();
        for m in removals {
            shard.remove(&m);
        }
    }
}

fn log_mem(step: &str) {
    if let Ok(statm) = std::fs::read_to_string("/proc/self/statm") {
        if let Some(rss_pages) = statm.split_whitespace().nth(1) {
            if let Ok(pages) = rss_pages.parse::<usize>() {
                let rss_mb = (pages * 4096) as f64 / 1_048_576.0;
                eprintln!("[MEM] {:.2} MB | {}", rss_mb, step);
            }
        }
    }
}