use std::collections::HashMap;
use std::io::{self, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::fs::File;
use anyhow::Result;
use rayon::prelude::*;
use crate::database::{AccessionRegistry, KmerDatabase, Hit};
use crate::minimizer::MinimizerScanner;
use rustc_hash::FxHashMap;

pub struct QueryConfig {
    pub db_prefix: String,
    pub read1_file: String,
    pub read2_file: Option<String>,
    pub threads: usize,
    pub use_accessions: bool,
    pub coverage_threshold: f64,
    pub output_prefix: String,
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

    let acc_path = format!("{}.accessions", config.db_prefix);
    let acc_registry = if config.use_accessions {
        match AccessionRegistry::load(&acc_path) {
            Ok(reg) => Some(reg),
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

    let is_paired = config.read2_file.is_some();

    eprintln!("Classifying reads");
    if is_paired {
        eprintln!("- Detected paired-end reads");
    }
    let start = std::time::Instant::now();

    let classified = AtomicUsize::new(0);
    let unclassified = AtomicUsize::new(0);

    let records1 = read_sequences(&config.read1_file)?;
    let records2 = if let Some(ref r2) = config.read2_file {
        let r2_recs = read_sequences(r2)?;
        if r2_recs.len() != records1.len() {
            anyhow::bail!(
                "Read count mismatch: R1 has {} reads, R2 has {} reads",
                records1.len(),
                r2_recs.len()
            );
        }
        Some(r2_recs)
    } else {
        None
    };

    let output_chunks: Mutex<Vec<Vec<u8>>> = Mutex::new(Vec::new());
    let report_counts: Mutex<HashMap<u32, usize>> = Mutex::new(HashMap::new());
    let batch_size = 5000;

    let num_records = records1.len();
    let batch_ranges: Vec<(usize, usize)> = (0..num_records)
        .step_by(batch_size)
        .map(|start| (start, (start + batch_size).min(num_records)))
        .collect();

    batch_ranges.par_iter().for_each(|&(batch_start, batch_end)| {
        let mut local_output: Vec<u8> = Vec::with_capacity(batch_size * 200);
        let mut local_report: HashMap<u32, usize> = HashMap::new();

        // Reusable buffers — allocated once per batch, reused every read
        let mut hits1: Vec<Hit> = Vec::new();
        let mut hits2: Vec<Hit> = Vec::new();
        let mut minimizer_buf1: Vec<u64> = Vec::new();
        let mut minimizer_buf2: Vec<u64> = Vec::new();
        let mut acc_counts: FxHashMap<u32, usize> = FxHashMap::default();

        for i in batch_start..batch_end {
            let record = &records1[i];
            let seq1 = record.sequence.as_bytes();

            // Query R1 — clears and fills hits1
            hits1.clear();
            db.query_into(&scanner, seq1, &mut minimizer_buf1, &mut hits1);

            let (seq2_len, has_r2) = if let Some(ref r2) = records2 {
                let r2_rec = &r2[i];
                let seq2 = r2_rec.sequence.as_bytes();
                hits2.clear();
                db.query_into(&scanner, seq2, &mut minimizer_buf2, &mut hits2);
                (seq2.len(), true)
            } else {
                (0, false)
            };

            // Coverage calc
            let read1_window = if seq1.len() >= db.k() { seq1.len() - db.k() + 1 } else { 1 };
            let read2_window = if has_r2 {
                if seq2_len >= db.k() { seq2_len - db.k() + 1 } else { 1 }
            } else {
                0
            };
            let total_window = read1_window + read2_window;

            let mut valid_hits: usize = 0;
            let mut taxid_counts = [0u32; 256];
            acc_counts.clear();

            // Count hits from R1
            for hit in hits1.iter() {
                if hit.is_hit {
                    valid_hits += 1;
                    taxid_counts[hit.taxid_idx as usize] += 1;
                    if acc_registry.is_some() {
                        for &acc in hit.accessions {
                            *acc_counts.entry(acc).or_default() += 1;
                        }
                    }
                }
            }

            // Count hits from R2
            if has_r2 {
                for hit in hits2.iter() {
                    if hit.is_hit {
                        valid_hits += 1;
                        taxid_counts[hit.taxid_idx as usize] += 1;
                        if acc_registry.is_some() {
                            for &acc in hit.accessions {
                                *acc_counts.entry(acc).or_default() += 1;
                            }
                        }
                    }
                }
            }

            let coverage = valid_hits as f64 / total_window as f64;

            if valid_hits == 0 || coverage < coverage_threshold {
                if is_paired {
                    let _ = write!(
                        local_output,
                        "U\t{}\t0\t{}|{}\t0:0\n",
                        record.header,
                        seq1.len(),
                        seq2_len,
                    );
                } else {
                    let _ = write!(
                        local_output,
                        "U\t{}\t0\t{}\t0:0\n",
                        record.header,
                        seq1.len()
                    );
                }
                unclassified.fetch_add(1, Ordering::Relaxed);
            } else {
                // Find best taxid from the fixed-size array
                let mut best_taxid_idx: u8 = 0;
                let mut best_count: u32 = 0;
                for (idx, &count) in taxid_counts.iter().enumerate() {
                    if count > best_count {
                        best_count = count;
                        best_taxid_idx = idx as u8;
                    }
                }

                let best_taxid = db.true_taxid(best_taxid_idx);

                if is_paired {
                    let _ = write!(
                        local_output,
                        "C\t{}\t{}\t{}|{}\t",
                        record.header, best_taxid, seq1.len(), seq2_len
                    );
                } else {
                    let _ = write!(
                        local_output,
                        "C\t{}\t{}\t{}\t",
                        record.header, best_taxid, seq1.len()
                    );
                }

                // Write R1 hit pattern
                write_hit_pattern(&mut local_output, &hits1, &db);

                if has_r2 {
                    local_output.extend_from_slice(b" |:| ");
                    write_hit_pattern(&mut local_output, &hits2, &db);
                }

                // Accession info
                if acc_registry.is_some() && !acc_counts.is_empty() {
                    local_output.push(b'\t');
                    let reg = acc_registry.as_ref().unwrap();
                    for (acc_id, count) in acc_counts.iter() {
                        let _ = write!(
                            local_output,
                            "{}:{} ",
                            reg.get_name(*acc_id),
                            count
                        );
                    }
                }

                let _ = write!(local_output, "\tcov:{:.3}\n", coverage);
                *local_report.entry(best_taxid).or_default() += 1;
                classified.fetch_add(1, Ordering::Relaxed);
            }
        }

        let mut chunks = output_chunks.lock().unwrap();
        chunks.push(local_output);

        let mut global_report = report_counts.lock().unwrap();
        for (taxid, count) in local_report {
            *global_report.entry(taxid).or_default() += count;
        }
    });

    // Write output
    {
        let f = File::create(format!("{}.output", config.output_prefix))?;
        let mut writer = io::BufWriter::new(f);
        let chunks = output_chunks.into_inner().unwrap();
        for chunk in chunks {
            writer.write_all(&chunk)?;
        }
        writer.flush()?;
    }

    // Write report
    {
        let f = File::create(format!("{}.report", config.output_prefix))?;
        let mut writer = io::BufWriter::new(f);

        let counts = report_counts.into_inner().unwrap();
        let mut sorted: Vec<(u32, usize)> = counts.into_iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(&a.1));

        writeln!(writer, "Target_TaxID\tRead_Count\tRatio")?;
        for (taxid, count) in &sorted {
            let ratio = if num_records as f64 > 0.0 { *count as f64 / num_records as f64 } else { 0.0 };
            writeln!(writer, "{}\t{}\t{:0.3}", taxid, count, ratio)?;
        }
        writer.flush()?;
    }

    let elapsed = start.elapsed();
    let c = classified.load(Ordering::Relaxed);
    let u = unclassified.load(Ordering::Relaxed);

    if c + u != num_records {
        eprintln!("Warning C+U is not adding up to number of reads");
    }

    eprintln!("Classified:   {}", c);
    eprintln!("Unclassified: {}", u);
    eprintln!("Total:        {}", num_records);
    eprintln!("Time:         {:.2}s", elapsed.as_secs_f64());
    if elapsed.as_secs_f64() > 0.0 {
        eprintln!(
            "Throughput:   {:.0} reads/s",
            num_records as f64 / elapsed.as_secs_f64()
        );
    }

    Ok(())
}

/// Write run-length encoded hit pattern for one end directly into a byte buffer.
fn write_hit_pattern(out: &mut Vec<u8>, hits: &[Hit], db: &KmerDatabase) {
    if hits.is_empty() {
        out.extend_from_slice(b"0:0");
        return;
    }

    let mut current_taxid = if hits[0].is_hit {
        db.true_taxid(hits[0].taxid_idx)
    } else {
        0
    };
    let mut run_len: u32 = 1;

    for hit in hits.iter().skip(1) {
        let t = if hit.is_hit {
            db.true_taxid(hit.taxid_idx)
        } else {
            0
        };
        if t == current_taxid {
            run_len += 1;
        } else {
            let _ = write!(out, "{}:{} ", current_taxid, run_len);
            current_taxid = t;
            run_len = 1;
        }
    }
    let _ = write!(out, "{}:{}", current_taxid, run_len);
}

fn read_sequences(path: &str) -> Result<Vec<FastqRecord>> {
    use needletail::parse_fastx_file;

    let mut records = Vec::new();
    let mut reader = parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Cannot open reads file {}", e))?;

    while let Some(result) = reader.next() {
        let rec = result.map_err(|e| anyhow::anyhow!("Not sure what is going on, but cant read records in {}", e))?;

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