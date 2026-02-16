use std::collections::HashMap;
use std::io::{self, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use anyhow::Result;
use rayon::prelude::*;

use crate::database::{AccessionRegistry, KmerDatabase};
use crate::minimizer::MinimizerScanner;

pub struct QueryConfig {
    pub db_prefix: String,
    pub reads_file: String,
    pub threads: usize,
    pub use_accessions: bool,
    pub is_paired: bool,
    pub coverage_threshold: f64,
}

struct FastqRecord {
    header: String,
    sequence: String,
}

pub fn run_query(config: QueryConfig) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok();
    let db = KmerDatabase::load(&config.db_prefix)?;
    let coverage_threshold = config.coverage_threshold;

    // Load accession registry if available and requested
    let acc_path = format!("{}.accessions", config.db_prefix);
    let acc_registry = if config.use_accessions {
        match AccessionRegistry::load(&acc_path) {
            Ok(reg) => {
                Some(reg)
            }
            Err(_) => {
                eprintln!("No accession registry found, skipping accession output");
                None
            }
        }
    } else {
        None
    };

    let scanner = MinimizerScanner::new(
        db.k(),
        db.l(),
        db.spaced_seed_mask,
        db.toggle_mask,
    );

    eprintln!("Classifying reads from {}...", config.reads_file);
    let start = std::time::Instant::now();

    let classified = AtomicUsize::new(0);
    let unclassified = AtomicUsize::new(0);

    // use needletail
    let records = read_sequences(&config.reads_file)?;

    let output_chunks: Mutex<Vec<Vec<u8>>> = Mutex::new(Vec::new());

    let batch_size = 5000;
    records.par_chunks(batch_size).for_each(|batch| {
        let mut local_output = String::with_capacity(batch_size * 200);

        for record in batch {
            let seq = record.sequence.as_bytes();
            let hits = db.query(&scanner, seq);

            let max_k = if seq.len() >= db.k() {
                seq.len() - db.k() + 1
            } else {
                1
            };

            let mut valid_hits = 0;
            let mut taxid_counts: HashMap<u8, usize> = HashMap::new();
            let mut acc_counts: HashMap<u32, usize> = HashMap::new();

            for hit in &hits {
                if hit.is_hit {
                    valid_hits += 1;
                    *taxid_counts.entry(hit.taxid_idx).or_default() += 1;
                    if acc_registry.is_some() {
                        for &acc in &hit.accessions {
                            *acc_counts.entry(acc).or_default() += 1;
                        }
                    }
                }
            }

            let coverage = valid_hits as f64 / max_k as f64;

            if valid_hits == 0 || coverage < coverage_threshold {
                local_output.push_str(&format!(
                    "U\t{}\t0\t{}\t0:0\n",
                    record.header,
                    seq.len()
                ));
                unclassified.fetch_add(1, Ordering::Relaxed);
            } else {
                
                

                // Find best taxid
                let best_taxid_idx = taxid_counts
                    .iter()
                    .max_by_key(|(_, &count)| count)
                    .map(|(&idx, _)| idx)
                    .unwrap_or(0);

                let best_taxid = db.true_taxid(best_taxid_idx);

                // Build classification string
                local_output.push_str(&format!(
                    "C\t{}\t{}\t{}\t",
                    record.header, best_taxid, seq.len()
                ));

                // Run-length encode the hit pattern
                if !hits.is_empty() {
                    let mut current_taxid = if hits[0].is_hit {
                        db.true_taxid(hits[0].taxid_idx)
                    } else {
                        0
                    };
                    let mut run_len = 1;

                    for hit in hits.iter().skip(1) {
                        let t = if hit.is_hit {
                            db.true_taxid(hit.taxid_idx)
                        } else {
                            0
                        };
                        if t == current_taxid {
                            run_len += 1;
                        } else {
                            local_output.push_str(&format!("{}:{} ", current_taxid, run_len));
                            current_taxid = t;
                            run_len = 1;
                        }
                    }
                    local_output.push_str(&format!("{}:{}", current_taxid, run_len));
                }

                // Accession info
                if acc_registry.is_some() && !acc_counts.is_empty() {
                    local_output.push('\t');
                    let reg = acc_registry.as_ref().unwrap();
                    for (acc_id, count) in &acc_counts {
                        local_output.push_str(&format!(
                            "{}:{} ",
                            reg.get_name(*acc_id),
                            count
                        ));
                    }
                }

                // Coverage
                local_output.push_str(&format!("\tcov:{:.3}\n", coverage));

                classified.fetch_add(1, Ordering::Relaxed);
            }
        }

        // Collect batch output
        let mut chunks = output_chunks.lock().unwrap();
        chunks.push(local_output.into_bytes());
    });

    // Write all output
    {
        let mut stdout = io::BufWriter::new(io::stdout().lock());
        let chunks = output_chunks.into_inner().unwrap();
        for chunk in chunks {
            stdout.write_all(&chunk)?;
        }
        stdout.flush()?;
    }

    let elapsed = start.elapsed();
    let c = classified.load(Ordering::Relaxed);
    let u = unclassified.load(Ordering::Relaxed);
    let total = c + u;

    eprintln!("\nResult");
    eprintln!("Classified:   {}", c);
    eprintln!("Unclassified: {}", u);
    eprintln!("Total:        {}", total);
    eprintln!(
        "Rate:         {:.1}%",
        if total > 0 {
            100.0 * c as f64 / total as f64
        } else {
            0.0
        }
    );
    eprintln!("Time:         {:.2}s", elapsed.as_secs_f64());
    if elapsed.as_secs_f64() > 0.0 {
        eprintln!(
            "Throughput:   {:.0} reads/s",
            total as f64 / elapsed.as_secs_f64()
        );
    }

    Ok(())
}

fn read_sequences(path: &str) -> Result<Vec<FastqRecord>> {
    use needletail::parse_fastx_file;

    let mut records = Vec::new();
    let mut reader = parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Cannot open reads file: {}", e))?;

    while let Some(result) = reader.next() {
        let rec = result.map_err(|e| anyhow::anyhow!("Error reading record: {}", e))?;

        let header = std::str::from_utf8(rec.id())
            .unwrap_or("")
            .to_string();

        let sequence = std::str::from_utf8(&rec.seq())
            .unwrap_or("")
            .to_uppercase();

        if !sequence.is_empty() {
            records.push(FastqRecord { header, sequence });
        }
    }

    Ok(records)
}