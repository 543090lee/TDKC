use std::collections::{BTreeSet, HashMap, HashSet};
use std::sync::{Arc, RwLock};
use std::time::Instant;
use anyhow::Result;
use boomphf::Mphf;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use crate::database::{AccessionRegistry, EqClassAccessions, KmerDatabase};
use crate::minimizer::{create_spaced_seed_mask, MinimizerScanner, SPACED_PATTERN, TOGGLE_MASK};
use crate::prep::load_prelim_map;
use crate::taxonomy::{load_target_taxids, BfsTaxonomy, TargetTaxIDManager, TaxonomyTree};
use crate::utils::{init_thread_pool,segment_ranges};
use crate::hash::compute_fingerprint;

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
    pub names_dmp_path: String,
}

const NUM_SHARDS: usize = 256;
const SHARD_SHIFT: usize = 56;
const BATCH_ITEM_LIMIT: usize = 2048;
const BATCH_BYTE_LIMIT: usize = 256 * 1024 * 1024;
const PAR_CHUNK_SIZE: usize = 64;
const SEG_TARGET_LEN: usize = 4 * 1024 * 1024;

pub fn run_build(config: BuildConfig) -> Result<()> {
    let spaced_seed_mask = create_spaced_seed_mask(SPACED_PATTERN);
    init_thread_pool(config.threads);

    let total_start = Instant::now();

    eprintln!("Loading taxonomy...");
    let taxonomy = TaxonomyTree::load(&config.nodes_file)?;

    eprintln!("Loading target TaxIDs...");
    let targets = load_target_taxids(&config.targets_file)?;
    let taxid_manager = TargetTaxIDManager::new(&targets, &taxonomy);

    // Roll up 
    let relevant_taxids = taxid_manager.all_relevant_taxids();
    // intenral taxIDs
    let bfs_tax = BfsTaxonomy::build(&taxonomy, &relevant_taxids);
    // taxID to name matcher
    let name_file = format!("{}.names", config.db_prefix);
    crate::taxonomy::lookup_names(&config.names_dmp_path, &targets, &name_file)?;
    
    eprintln!("\nLoading prelim_map...");
    let acc_to_taxid = load_prelim_map(&config.prelim_map_file)?;

    let acc_to_internal: HashMap<String, u32> = acc_to_taxid
        .iter()
        .filter_map(|(acc, ext_taxid): (&String, &u32)| {
            bfs_tax.to_internal(*ext_taxid).map(|int_id| (acc.clone(), int_id))
        })
        .collect();
    eprintln!(
        "{} accessions have taxids in taxonomy tree",
        acc_to_internal.len()
    );

    let target_accessions: HashSet<String> = acc_to_internal
        .iter()
        .filter(|(_, int_id): &(&String, &u32)| bfs_tax.is_relevant(**int_id))
        .map(|(acc, _): (&String, &u32)| acc.clone())
        .collect();
    
        eprintln!(
        "{} accessions belong to target clades",
        target_accessions.len()
    );

    let scanner = MinimizerScanner::new(config.k, config.l, spaced_seed_mask, TOGGLE_MASK);

    let mut accession_registry = if config.track_accessions {
        let reg = presort_accessions_by_lineage(
            &target_accessions,
            &acc_to_taxid,
            &taxonomy,
            &taxid_manager,
        )?;
        Some(reg)
    } else {
        None
    };

    eprintln!("\nExtracting target minimizers...");
    let global = ShardedMinimizerMap::new();

    extract_target_minimizers(
        &config.target_fasta_file,
        &target_accessions,
        &acc_to_internal,
        &bfs_tax,
        &scanner,
        &global,
        config.track_accessions,
        &mut accession_registry,
    )?;

    let total_target_minimizers: usize = global
        .shards
        .iter()
        .map(|s| s.read().unwrap().len())
        .sum();

    eprintln!("\nChallenging minimizers...");
    let to_remove = collect_background_hits(
        &config.fasta_file,
        &target_accessions,
        &scanner,
        &global,
    )?;

    let remove_count = to_remove.len();
    eprintln!("{} minimizers found in non-target, removing...", remove_count);

    challenge_bulk_remove(&global, &to_remove);
    drop(to_remove);

    let total_surviving: usize = global
        .shards
        .iter()
        .map(|s| s.read().unwrap().len())
        .sum();
    eprintln!(
        "{} minimizers survived, abd {} removed by challenge",
        total_surviving,
        total_target_minimizers.saturating_sub(total_surviving)
    );

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

    db.save(&config.db_prefix)?;

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

    eprintln!("Found {} unique taxIDs", unique_taxids.len());
    eprintln!("{} unique minimizers", num_minimizers);

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
    let mut keys: Vec<u64> = Vec::with_capacity(num_minimizers);
    for shard_lock in &global.shards {
        let shard = shard_lock.read().unwrap();
        for (&minimizer, entry) in shard.iter() {
            if bfs_tax.is_relevant(entry.int_taxid) {
                keys.push(minimizer);
            }
        }
    }

    eprintln!("Building MPHF...");
    let mphf = Mphf::new(2.0, &keys);
    drop(keys);

    eprintln!("\nPopulating arrays...");
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
    }

    drop(hash_to_class);

    if let Some(ref mut acc) = accessions {
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
        None,
    ))
}

fn presort_accessions_by_lineage(
    target_accessions: &HashSet<String>,
    acc_to_taxid: &FxHashMap<String, u32>,
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

fn extract_target_minimizers(
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

        while let Some(record) = reader.next() {
            let rec = record.map_err(|e| anyhow::anyhow!("FASTA error: {}", e))?;
            let accession = crate::utils::extract_accession(rec.id());

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

            let ranges = segment_ranges(seq_len, k, SEG_TARGET_LEN);

            for (start, end) in ranges {
                batch_bytes += end - start;
                work_items.push(WorkItem {
                    seq_data: Arc::clone(&seq_data),
                    start,
                    end,
                    int_taxid,
                });
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
        }

        if !work_items.is_empty() || !acc_records.is_empty() {
            let _ = tx.send((work_items, acc_records));
        }
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
        Err(_) => anyhow::bail!("Reader thread panicked"),
    }

    Ok(())
}
struct ChallengeSegment {
    seq_data: Arc<Vec<u8>>,
    start: usize,
    end: usize,
}

impl ChallengeSegment {
    #[inline]
    fn seq_slice(&self) -> &[u8] {
        &self.seq_data[self.start..self.end]
    }
}

fn collect_background_hits(
    fasta_path: &str,
    target_accessions: &HashSet<String>,
    scanner: &MinimizerScanner,
    global: &ShardedMinimizerMap,
) -> Result<Vec<u64>> {
    let k = scanner.k();

    let (tx, rx) = crossbeam_channel::bounded::<Vec<ChallengeSegment>>(4);
    let fasta_file = fasta_path.to_string();
    let target_acc = target_accessions.clone();
    let reader_handle = std::thread::spawn(move || -> Result<()> {
        let mut reader = needletail::parse_fastx_file(&fasta_file)
            .map_err(|e| anyhow::anyhow!("Cannot open FASTA: {}", e))?;

        let mut work_items: Vec<ChallengeSegment> = Vec::new();
        let mut batch_bytes = 0usize;

        while let Some(record) = reader.next() {
            let rec = record.map_err(|e| anyhow::anyhow!("FASTA error: {}", e))?;
            let accession = crate::utils::extract_accession(rec.id());

            if target_acc.contains(accession) {                
                continue;
            }

            let seq_data = Arc::new(rec.seq().into_owned());
            let seq_len = seq_data.len();
            let ranges = segment_ranges(seq_len, k, SEG_TARGET_LEN);
            for (start, end) in ranges {
                batch_bytes += end - start;
                work_items.push(ChallengeSegment {
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
        }
        if !work_items.is_empty() {
            let _ = tx.send(work_items);
        }
        Ok(())
    });

    let mut all_hits: Vec<u64> = Vec::new();

    for work_items in rx {
        let batch_hits: Vec<u64> = work_items
            .par_chunks(PAR_CHUNK_SIZE)
            .fold(
                || Vec::new(),
                |mut local_hits: Vec<u64>, chunk: &[ChallengeSegment]| {
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

fn challenge_bulk_remove(global: &ShardedMinimizerMap, to_remove: &[u64]) {
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
