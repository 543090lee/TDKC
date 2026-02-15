use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::Mutex;

use anyhow::Result;
use rayon::prelude::*;

use crate::database::{AccessionRegistry, ExtractedKmer, KmerDatabaseBuilder};
use crate::fasta_index::FastaIndex;
use crate::minimizer::create_spaced_seed_mask;
use crate::taxonomy::{load_target_taxids, TargetTaxIDManager, TaxonomyTree};

pub struct BuildConfig {
    pub kraken_file: String,
    pub fasta_file: String,
    pub targets_file: String,
    pub nodes_file: String,
    pub db_prefix: String,
    pub threads: usize,
    pub track_accessions: bool,
}

pub fn run_build(config: BuildConfig) -> Result<()> {
    const K: usize = 35;
    const L: usize = 31;
    let spaced_pattern = "1111111111111111111110101010101";
    let spaced_seed_mask = create_spaced_seed_mask(spaced_pattern);
    let toggle_mask: u64 = 0;

    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok();

    let total_start = std::time::Instant::now();

    eprintln!("Loading Taxonomy");
    let taxonomy = TaxonomyTree::load(&config.nodes_file)?;

    eprintln!("\nLoading Target TaxIDs");
    let targets = load_target_taxids(&config.targets_file)?;
    let taxid_manager = TargetTaxIDManager::new(&targets, &taxonomy);

    eprintln!("\nIndexing FASTA");
    let fasta_index = FastaIndex::new(&config.fasta_file)?;

    eprintln!("\nExtracting K-mers");
    let relevant_taxids = taxid_manager.all_relevant_taxids();

    let accession_registry = if config.track_accessions {
        Some(Mutex::new(AccessionRegistry::new()))
    } else {
        None
    };

    let kmers = extract_kmers_from_kraken(
        &config.kraken_file,
        &fasta_index,
        &taxid_manager,
        &relevant_taxids,
        &accession_registry,
        K,
    )?;

    eprintln!("Extracted {} unique k-mers total", kmers.len());

    eprintln!("\nBuilding Database");
    let builder = KmerDatabaseBuilder::new(K, L, spaced_seed_mask, toggle_mask, config.track_accessions);
    let db = builder.build(&kmers, config.threads)?;
    db.save(&config.db_prefix)?;

    // Save accession registry if tracking
    if let Some(ref reg_mutex) = accession_registry {
        let reg = reg_mutex.lock().unwrap();
        reg.save(&format!("{}.accessions", config.db_prefix))?;
        eprintln!("Saved {} accessions", reg.len());
    }

    let total_elapsed = total_start.elapsed();
    eprintln!("  Done!");
    eprintln!("Took {:.2}s", total_elapsed.as_secs_f64());

    Ok(())
}

fn extract_kmers_from_kraken(
    kraken_file: &str,
    fasta_index: &FastaIndex,
    taxid_manager: &TargetTaxIDManager,
    relevant_taxids: &HashSet<u32>,
    accession_registry: &Option<Mutex<AccessionRegistry>>,
    k: usize,
) -> Result<Vec<ExtractedKmer>> {
    let file = File::open(kraken_file)?;
    let reader = BufReader::with_capacity(1024 * 1024, file);

    // Read all lines into batches for parallel processing
    let lines: Vec<String> = reader.lines().collect::<std::io::Result<Vec<_>>>()?;
    eprintln!("Processing {} reads...", lines.len());

    let batch_size = 50_000;
    let all_kmers: Mutex<Vec<ExtractedKmer>> = Mutex::new(Vec::new());
    let global_seen: Mutex<HashSet<String>> = Mutex::new(HashSet::new());

    lines
        .par_chunks(batch_size)
        .for_each(|batch| {
            let mut local_kmers: Vec<ExtractedKmer> = Vec::new();
            let mut local_seen: HashSet<String> = HashSet::new();

            for line in batch {
                if let Some(mut extracted) = process_kraken_line(
                    line,
                    fasta_index,
                    taxid_manager,
                    relevant_taxids,
                    accession_registry,
                    k,
                    &mut local_seen,
                ) {
                    local_kmers.append(&mut extracted);
                }
            }

            // Merge into global
            let mut global = global_seen.lock().unwrap();
            let mut all = all_kmers.lock().unwrap();
            for kmer in local_kmers {
                if global.insert(kmer.sequence.clone()) {
                    all.push(kmer);
                }
            }
        });

    Ok(all_kmers.into_inner().unwrap())
}

fn process_kraken_line(
    line: &str,
    fasta_index: &FastaIndex,
    taxid_manager: &TargetTaxIDManager,
    relevant_taxids: &HashSet<u32>,
    accession_registry: &Option<Mutex<AccessionRegistry>>,
    k: usize,
    local_seen: &mut HashSet<String>,
) -> Option<Vec<ExtractedKmer>> {
    let mut parts = line.split('\t');
    let _classification = parts.next()?;
    let seq_id = parts.next()?;
    let _taxid_str = parts.next()?;
    let _seq_len = parts.next()?;
    let lca_mapping = parts.next().unwrap_or("");

    let sequence = fasta_index.get_sequence(seq_id)?;
    if sequence.is_empty() {
        return None;
    }

    let acc_id = accession_registry.as_ref().map(|reg| {
        reg.lock().unwrap().get_or_create(seq_id)
    });

    let seq_bytes = sequence.as_bytes();
    let seq_len = seq_bytes.len();
    let mut kmers = Vec::new();
    let mut kmer_read_index: usize = 0;

    for part in lca_mapping.split_whitespace() {
        let colon_pos = match part.find(':') {
            Some(p) => p,
            None => continue,
        };

        let taxid_part = &part[..colon_pos];
        let count_part = &part[colon_pos + 1..];

        if taxid_part == "cov" {
            continue;
        }

        let taxid: u32 = match taxid_part.parse() {
            Ok(v) => v,
            Err(_) => {
                if let Ok(count) = count_part.parse::<usize>() {
                    kmer_read_index += count;
                }
                continue;
            }
        };

        let count: usize = match count_part.parse() {
            Ok(v) => v,
            Err(_) => continue,
        };

        if relevant_taxids.contains(&taxid) {
            if let Some(target_taxid) = taxid_manager.get_target(taxid) {
                for j in 0..count {
                    let start = kmer_read_index + j;
                    let end = start + k;
                    if end <= seq_len {
                        let kmer_str = &sequence[start..end];

                        // Validate: only ACGT
                        let valid = kmer_str
                            .bytes()
                            .all(|c| matches!(c, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'));

                        if valid && !local_seen.contains(kmer_str) {
                            local_seen.insert(kmer_str.to_string());
                            kmers.push(ExtractedKmer {
                                sequence: kmer_str.to_string(),
                                taxid: target_taxid,
                                accession_id: acc_id,
                            });
                        }
                    }
                }
            }
        }

        kmer_read_index += count;
    }

    if kmers.is_empty() {
        None
    } else {
        Some(kmers)
    }
}