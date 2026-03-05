use std::collections::HashMap;
use std::io::{self, Write};
use std::fs::File;
use anyhow::Result;
use rayon::prelude::*;
use crate::database::{AccessionRegistry, KmerDatabase, Hit};
use crate::minimizer::MinimizerScanner;
use rustc_hash::FxHashMap;
use itertools::Itertools;

pub struct QueryConfig {
    pub db_prefix: String,
    pub read1_file: String,
    pub read2_file: Option<String>,
    pub threads: usize,
    pub use_accessions: bool,
    pub coverage_threshold: Option<f64>,
    pub output_prefix: String,
}

// Not using string anymore here, now to u8
struct FastqRecord {
    header: Vec<u8>,
    sequence: Vec<u8>,
}

struct ReadBatch {
    records1: Vec<FastqRecord>,
    records2: Option<Vec<FastqRecord>>,
}

// this goes to writer thread
struct BatchResult {
    output_data: Vec<u8>,
    report_counts: FxHashMap<u32, usize>,
    classified: usize,
    unclassified: usize,
    num_records: usize,
}

pub fn run_query(config: QueryConfig) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.threads)
        .build_global()
        .ok();

    let db_start = std::time::Instant::now();

    let db = KmerDatabase::load(&config.db_prefix, config.use_accessions)?;

    let acc_path = format!("{}.accessions", config.db_prefix);
    let acc_registry = if config.use_accessions {
        match AccessionRegistry::load_for_query(&acc_path) {
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

    let db_elapsed = db_start.elapsed();

    let is_paired = config.read2_file.is_some();

    if is_paired {
        eprintln!("Running paired-end reads");
    } else {
        eprintln!("Running single-end reads");
    }

    let classify_start = std::time::Instant::now();

    // thread only to write, no lock, so less overhead
    let output_file = File::create(format!("{}.output", config.output_prefix))?;

    // we need sync channels here, since rayon worker threads will be much faster and more, just spitting
    // BatchResults into writer thread, but if we dont have sync channels, rayon threads will just stall...

    let (writer_tx, writer_rx) = std::sync::mpsc::sync_channel::<BatchResult>(config.threads * 4);

    // literally writer thread will act as a mutex, but without contention
    let writer_handle = std::thread::spawn(move || -> Result<(HashMap<u32, usize>, usize, usize, usize)> {

        //now we dont have to worry about flushing it out, since BufWriter will automatically write to disk and reset
        //when it;s full
        let mut writer = io::BufWriter::with_capacity(4 * 1024 * 1024, output_file);
        let mut global_report: HashMap<u32, usize> = HashMap::new();

        // I could use AtomicUsize for these, since these are single values, but too lazy
        let mut total_c: usize = 0;
        let mut total_u: usize = 0;
        let mut total_n: usize = 0;

        for batch in writer_rx {
            writer.write_all(&batch.output_data)?;
            for (taxid, count) in batch.report_counts {
                *global_report.entry(taxid).or_default() += count;
            }
            total_c += batch.classified;
            total_u += batch.unclassified;
            total_n += batch.num_records;
        }
        writer.flush()?;
        Ok((global_report, total_c, total_u, total_n))
    });


    let batch_size = 5_000;

    // by doing this, reader keep pushes batches into crossbeam channel, and bridge lets each rayon worker
    // grab one batch whenever it's free (mpmc), so no contention like before, where if one thread slow, then all had 
    // to wait before starting on the nextc chunk (mpsc)

    let (batch_tx, batch_rx) = crossbeam_channel::bounded::<ReadBatch>(config.threads * 2);

    // Reader thread(s): produce batches directly (no intermediate chunks)
    let r1_path = config.read1_file.clone();
    let r2_path = config.read2_file.clone();

    let reader_handle = std::thread::spawn(move || -> Result<()> {
        if let Some(ref r2_path_str) = r2_path {
            read_paired_batches(&r1_path, r2_path_str, batch_size, &batch_tx)?;
        } else {
            read_single_batches(&r1_path, batch_size, &batch_tx)?;
        }
        // sender is done, so done reading (EOF)
        drop(batch_tx); 
        Ok(())
    });

    //par_bridge() let's channel to give batch to rayon thread whenever it's free
    batch_rx.into_iter().par_bridge().for_each(|batch| {
        let result = classify_batch(
            &batch,
            &db,
            &scanner,
            &acc_registry,
            &config.coverage_threshold,
            is_paired,
            db.k,
            db.l
        );
        let _ = writer_tx.send(result);
    });

    drop(writer_tx);
    let (global_report, c, u, num_records) = match writer_handle.join() {
        Ok(Ok(v)) => v,
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("Writer thread panicked"),
    };

    match reader_handle.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("Reader thread panicked"),
    }

    let classify_elapsed = classify_start.elapsed();

    // Write report
    {
        let f = File::create(format!("{}.report", config.output_prefix))?;
        let mut writer = io::BufWriter::new(f);

        let mut sorted: Vec<(u32, usize)> = global_report.into_iter().collect();
        sorted.sort_by(|a, b| b.1.cmp(&a.1));

        writeln!(writer, "Target_TaxID\tRead_Count\tRatio")?;
        for (taxid, count) in &sorted {
            let ratio = if num_records > 0 { *count as f64 / num_records as f64 } else { 0.0 };
            writeln!(writer, "{}\t{}\t{:0.3}", taxid, count, ratio)?;
        }
        writer.flush()?;
    }

    if c + u != num_records {
        eprintln!("Warning C+U is not adding up to number of reads");
    }

    eprintln!("Classified:   {}", c);
    eprintln!("Unclassified: {}", u);
    eprintln!("Total:        {}", num_records);
    eprintln!("DB load:      {:.2}s", db_elapsed.as_secs_f64());
    eprintln!("Classify:     {:.2}s", classify_elapsed.as_secs_f64());
    eprintln!("Total time:   {:.2}s", db_elapsed.as_secs_f64() + classify_elapsed.as_secs_f64());
    if classify_elapsed.as_secs_f64() > 0.0 {
        eprintln!(
            "Throughput:   {:.0} reads/s",
            num_records as f64 / classify_elapsed.as_secs_f64()
        );
    }
    Ok(())
}

// Tried not to use mutex here since it was keep giving me extra lock contention overhead
// This was the main reason I switched to having a main writer thread

#[inline(never)]
fn classify_batch(
    batch: &ReadBatch,
    db: &KmerDatabase,
    scanner: &MinimizerScanner,
    acc_registry: &Option<AccessionRegistry>,
    coverage_threshold: &Option<f64>,
    is_paired: bool,
    k: usize,
    l: usize,

) -> BatchResult {
    let batch_len = batch.records1.len();
    let mut local_output: Vec<u8> = Vec::with_capacity(batch_len * 200);
    let mut local_report: FxHashMap<u32, usize> = FxHashMap::default();
    let mut local_classified: usize = 0;
    let mut local_unclassified: usize = 0;

    let mut hits1: Vec<Hit> = Vec::new();
    let mut hits2: Vec<Hit> = Vec::new();
    let mut minimizer_buf1: Vec<u64> = Vec::new();
    let mut minimizer_buf2: Vec<u64> = Vec::new();

    for i in 0..batch_len {
        let record = &batch.records1[i];
        let seq1 = &record.sequence;

        hits1.clear();
        db.query_into(scanner, seq1, &mut minimizer_buf1, &mut hits1);

        let (seq2_len, has_r2) = if let Some(ref r2) = batch.records2 {
            let r2_rec = &r2[i];
            let seq2 = &r2_rec.sequence;
            hits2.clear();
            db.query_into(scanner, seq2, &mut minimizer_buf2, &mut hits2);
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

        for hit in hits1.iter() {
            if hit.is_hit {
                valid_hits += 1;
                taxid_counts[hit.taxid_idx as usize] += 1;
            }
        }

        if has_r2 {
            for hit in hits2.iter() {
                if hit.is_hit {
                    valid_hits += 1;
                    taxid_counts[hit.taxid_idx as usize] += 1;
                }
            }
        }

        let coverage = valid_hits as f64 / total_window as f64;
        // I should add another threshold of minimu hit group that looks at different number of minimizer hits
        let threshold = coverage_threshold.unwrap_or((k as f64 - l as f64 + 1.0)*2.0 / (total_window as f64));
        if valid_hits == 0 || coverage <= threshold {
            if is_paired {
                let _ = write!(
                    local_output,
                    "U\t{}\t0\t{}|{}\t0:0\n",
                    unsafe { std::str::from_utf8_unchecked(&record.header) },
                    seq1.len(), seq2_len,
                );
            } else {
                let _ = write!(
                    local_output,
                    "U\t{}\t0\t{}\t0:0\n",
                    unsafe { std::str::from_utf8_unchecked(&record.header) },
                    seq1.len()
                );
            }
            local_unclassified += 1;
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
                    unsafe { std::str::from_utf8_unchecked(&record.header) },
                    best_taxid, seq1.len(), seq2_len
                );
            } else {
                let _ = write!(
                    local_output,
                    "C\t{}\t{}\t{}\t",
                    unsafe { std::str::from_utf8_unchecked(&record.header) },
                    best_taxid, seq1.len()
                );
            }

            write_hit_pattern(&mut local_output, &hits1, db, &acc_registry);

            if has_r2 {
                local_output.extend_from_slice(b" |:| ");
                write_hit_pattern(&mut local_output, &hits2, db, &acc_registry);
            }

            let _ = write!(local_output, "\tcov:{:.3}\n", coverage);
            *local_report.entry(best_taxid).or_default() += 1;
            local_classified += 1;
        }
    }

    BatchResult {
        output_data: local_output,
        report_counts: local_report,
        classified: local_classified,
        unclassified: local_unclassified,
        num_records: batch_len,
    }
}

// resolve accessions lazily from class_id only when writing output
fn write_hit_pattern(out: &mut Vec<u8>, hits: &[Hit], db: &KmerDatabase, acc_registry: &Option<AccessionRegistry>) {
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

    let mut current_class_id = hits[0].accession_class_id;
    let reg = acc_registry.as_ref();

    for hit in hits.iter().skip(1) {
        let t = if hit.is_hit {
            db.true_taxid(hit.taxid_idx)
        } else {
            0
        };
        if t == current_taxid {
            run_len += 1;
        } else {
            emit_run(out, current_taxid, run_len, current_class_id, reg, db);
            current_taxid = t;
            current_class_id = hit.accession_class_id;
            run_len = 1;
        }
    }
    emit_run(out, current_taxid, run_len, current_class_id, reg, db);
}

#[inline]
fn emit_run(
    out: &mut Vec<u8>,
    taxid: u32,
    run_len: u32,
    class_id: Option<u32>,
    reg: Option<&AccessionRegistry>,
    db: &KmerDatabase,
) {
    if let Some(r) = reg {
        if taxid != 0 {
            if let Some(cid) = class_id {
                // Only decode accessions here — this is the only place we need the actual names
                let accs = db.resolve_accessions(cid);
                if !accs.is_empty() {
                    let acc_str = accs.iter().map(|id| r.get_name(*id)).join(",");
                    let _ = write!(out, "{}-{{{}}}:{} ", taxid, acc_str, run_len);
                    return;
                }
            }
        }
    }
    let _ = write!(out, "{}:{} ", taxid, run_len);
}

fn read_single_batches(
    r1_path: &str,
    batch_size: usize,
    tx: &crossbeam_channel::Sender<ReadBatch>,
) -> Result<()> {
    let mut r1_reader = open_fastx(r1_path)?;
    loop {
        let recs1 = read_records(&mut r1_reader, batch_size)?;
        if recs1.is_empty() {
            break;
        }
        if tx.send(ReadBatch { records1: recs1, records2: None }).is_err() {
            break;
        }
    }
    Ok(())
}

fn read_paired_batches(
    r1_path: &str,
    r2_path: &str,
    batch_size: usize,
    tx: &crossbeam_channel::Sender<ReadBatch>,
) -> Result<()> {
    // R2 decompresses in its own thread, sends batches to us
    let (r2_tx, r2_rx) = crossbeam_channel::bounded::<Option<Vec<FastqRecord>>>(4);
    let r2_path_owned = r2_path.to_string();
    let bs = batch_size;

    let r2_handle = std::thread::spawn(move || -> Result<()> {
        let mut r2_reader = open_fastx(&r2_path_owned)?;
        loop {
            let recs2 = read_records(&mut r2_reader, bs)?;
            if recs2.is_empty() {
                let _ = r2_tx.send(None);
                break;
            }
            if r2_tx.send(Some(recs2)).is_err() {
                break;
            }
        }
        Ok(())
    });

    let mut r1_reader = open_fastx(r1_path)?;
    loop {
        let recs1 = read_records(&mut r1_reader, batch_size)?;
        if recs1.is_empty() {
            break;
        }
        let recs2 = match r2_rx.recv() {
            Ok(Some(r2)) => r2,
            Ok(None) => anyhow::bail!("R2 file ended before R1"),
            Err(_) => anyhow::bail!("R2 reader thread disconnected"),
        };
        if recs2.len() != recs1.len() {
            anyhow::bail!(
                "There is a read count mismatch: R1={}, R2={}",
                recs1.len(), recs2.len()
            );
        }
        if tx.send(ReadBatch { records1: recs1, records2: Some(recs2) }).is_err() {
            break;
        }
    }

    match r2_handle.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("R2 reader thread panicked"),
    }
    Ok(())
}

type FastxReader = Box<dyn needletail::FastxReader + Send>;

fn open_fastx(path: &str) -> Result<FastxReader> {
    let reader = needletail::parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Cannot open reads file {}", e))?;
    Ok(reader)
}

/// Read up to n records. Raw bytes — no to_uppercase(), no UTF-8 validation on sequence.
fn read_records(reader: &mut FastxReader, n: usize) -> Result<Vec<FastqRecord>> {
    let mut records = Vec::with_capacity(n);
    for _ in 0..n {
        match reader.next() {
            Some(Ok(rec)) => {
                let seq = rec.seq();
                if !seq.is_empty() {
                    records.push(FastqRecord {
                        header: rec.id().to_vec(),
                        sequence: seq.to_vec(),
                    });
                }
            }
            Some(Err(e)) => return Err(anyhow::anyhow!("Not sure what the error is but can't read records {}", e)),
            None => break,
        }
    }
    Ok(records)
}