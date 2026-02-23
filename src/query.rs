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

// A chunk of paired reads ready for classification
struct ReadChunk {
    records1: Vec<FastqRecord>,
    records2: Option<Vec<FastqRecord>>,
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

    if is_paired {
        eprintln!("Running paired-end reads");
    } else {
        eprintln!("Running single-end reads");
    }

    let start = std::time::Instant::now();

    let classified = AtomicUsize::new(0);
    let unclassified = AtomicUsize::new(0);
    let total_records = AtomicUsize::new(0);

    // Output file, batches write directly
    let output_file = File::create(format!("{}.output", config.output_prefix))?;
    let output_writer = Mutex::new(io::BufWriter::with_capacity(4 * 1024 * 1024, output_file));
    let report_counts: Mutex<HashMap<u32, usize>> = Mutex::new(HashMap::new());

    // we have double buffer here, where a reader thread keeps on reading
    // and the main thread gets the chunk and dispatch them into rayon for classification
    // 2 channels seem reasonable, since on 45M reads, rayon was little falling behind
    let chunk_size = 100_000;

    let (tx, rx) = std::sync::mpsc::sync_channel::<ReadChunk>(2);

    let r1_path = config.read1_file.clone();
    let r2_path = config.read2_file.clone();

    // Reader thread, not a rayon thread pool
    let reader_handle = std::thread::spawn(move || -> Result<()> {
        let mut r1_reader = open_fastx(&r1_path)?;
        let mut r2_reader = match r2_path {
            Some(ref p) => Some(open_fastx(p)?),
            None => None,
        };

        loop {
            let recs1 = read_chunk(&mut r1_reader, chunk_size)?;
            if recs1.is_empty() {
                break;
            }

            let recs2 = if let Some(ref mut r2) = r2_reader {
                let c2 = read_chunk(r2, chunk_size)?;
                if c2.len() != recs1.len() {
                    anyhow::bail!(
                        "Read count mismatch in chunk: R1 has {} reads, R2 has {} reads",
                        recs1.len(),
                        c2.len()
                    );
                }
                Some(c2)
            } else {
                None
            };

            if tx.send(ReadChunk { records1: recs1, records2: recs2 }).is_err() {
                break; // Receiver dropped
            }
        }
        Ok(())
    });

    // Classify chunks as they arrive 
    // Since the chunk size is 100K, if you do 10k batch size, there will only be 10 tasks
    // If you have more threads than that, maybe it's a good idea to go below. But going below that
    // plateus. At least locally on my laptop M4 Max
  
    let batch_size = 10000;

    for chunk in rx {
        let chunk_len = chunk.records1.len();
        total_records.fetch_add(chunk_len, Ordering::Relaxed);

        let batch_ranges: Vec<(usize, usize)> = (0..chunk_len)
            .step_by(batch_size)
            .map(|s| (s, (s + batch_size).min(chunk_len)))
            .collect();

        batch_ranges.par_iter().for_each(|&(batch_start, batch_end)| {
            let mut local_output: Vec<u8> = Vec::with_capacity(batch_size * 200);
            let mut local_report: HashMap<u32, usize> = HashMap::new();

            let mut hits1: Vec<Hit> = Vec::new();
            let mut hits2: Vec<Hit> = Vec::new();
            let mut minimizer_buf1: Vec<u64> = Vec::new();
            let mut minimizer_buf2: Vec<u64> = Vec::new();
            let mut acc_counts: FxHashMap<u32, usize> = FxHashMap::default();

            for i in batch_start..batch_end {
                let record = &chunk.records1[i];
                let seq1 = record.sequence.as_bytes();

                hits1.clear();
                db.query_into(&scanner, seq1, &mut minimizer_buf1, &mut hits1);

                let (seq2_len, has_r2) = if let Some(ref r2) = chunk.records2 {
                    let r2_rec = &r2[i];
                    let seq2 = r2_rec.sequence.as_bytes();
                    hits2.clear();
                    db.query_into(&scanner, seq2, &mut minimizer_buf2, &mut hits2);
                    (seq2.len(), true)
                } else {
                    (0, false)
                };

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
                            record.header, seq1.len(), seq2_len,
                        );
                    } else {
                        let _ = write!(
                            local_output,
                            "U\t{}\t0\t{}\t0:0\n",
                            record.header, seq1.len()
                        );
                    }
                    unclassified.fetch_add(1, Ordering::Relaxed);
                } else {
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

                    write_hit_pattern(&mut local_output, &hits1, &db);

                    if has_r2 {
                        local_output.extend_from_slice(b" |:| ");
                        write_hit_pattern(&mut local_output, &hits2, &db);
                    }

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

            // Write output directly
            {
                let mut writer = output_writer.lock().unwrap();
                let _ = writer.write_all(&local_output);
            }

            {
                let mut global_report = report_counts.lock().unwrap();
                for (taxid, count) in local_report {
                    *global_report.entry(taxid).or_default() += count;
                }
            }
        });
    }

    // Wait for reader thread and propagate errors
    match reader_handle.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("Reader thread panicked"),
    }

    // Flush output
    {
        let mut writer = output_writer.lock().unwrap();
        writer.flush()?;
    }

    let num_records = total_records.load(Ordering::Relaxed);

    // Write report
    {
        let f = File::create(format!("{}.report", config.output_prefix))?;
        let mut writer = io::BufWriter::new(f);

        let counts = report_counts.into_inner().unwrap();
        let mut sorted: Vec<(u32, usize)> = counts.into_iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(&a.1));

        writeln!(writer, "Target_TaxID\tRead_Count\tRatio")?;
        for (taxid, count) in &sorted {
            let ratio = if num_records > 0 { *count as f64 / num_records as f64 } else { 0.0 };
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

type FastxReader = Box<dyn needletail::FastxReader + Send>;

fn open_fastx(path: &str) -> Result<FastxReader> {
    let reader = needletail::parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Cannot open reads file {}", e))?;
    Ok(reader)
}

// Read up to n records from the reader, and returns empty vec at EOF
fn read_chunk(reader: &mut FastxReader, n: usize) -> Result<Vec<FastqRecord>> {
    let mut records = Vec::with_capacity(n);

    for _ in 0..n {
        match reader.next() {
            Some(Ok(rec)) => {
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
            Some(Err(e)) => return Err(anyhow::anyhow!("Not sure what the error is but can't read records {}", e)),
            None => break,
        }
    }

    Ok(records)
}