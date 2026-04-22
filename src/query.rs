use std::collections::HashMap;
use std::io::{self, Write};
use std::fs::File;
use anyhow::Result;
use rayon::prelude::*;
use crate::database::{AccessionRegistry, DomainBloomFilter, KmerDatabase, Hit};
use crate::minimizer::MinimizerScanner;
use rustc_hash::FxHashMap;
use itertools::Itertools;
use crate::utils::init_thread_pool;
use crate::taxonomy::{TAXID_ARCHAEA, TAXID_BACTERIA, TAXID_FUNGI, TAXID_VIRAL, resolve_domain_lca};
use crate::read_finder::Sample;

pub struct QueryConfig {
    pub db_prefix: String,
    pub samples: Vec<Sample>,
    pub threads: usize,
    pub use_accessions: bool,
    pub background: bool,
    pub minimum_hit_groups: usize,
    pub output_prefix: String,
    pub no_output: bool,
}

pub struct FlatBatch {
    pub data: Vec<u8>,
    pub offsets: Vec<(usize, usize, usize, usize)>,
}

pub struct ReadBatch {
    pub batch1: FlatBatch,
    pub batch2: Option<FlatBatch>,
}

// this goes to writer thread
struct BatchResult {
    output_data: Vec<u8>,
    report_counts: FxHashMap<u32, usize>,
    classified: usize,
    unclassified: usize,
    ambiguous: usize,
    num_records: usize,
}

pub fn run_query(config: QueryConfig) -> Result<()> {
    init_thread_pool(config.threads);

    let db = KmerDatabase::load(&config.db_prefix, config.use_accessions)?;

    let mut bg_filters = Vec::new();
    if config.background {
        let domains = [
            ("bacteria.bloom", TAXID_BACTERIA),
            ("archaea.bloom", TAXID_ARCHAEA),
            ("viral.bloom", TAXID_VIRAL),
            ("fungi.bloom", TAXID_FUNGI)
        ];
        for (ext, taxid) in domains {
            let path = format!("{}.{}", config.db_prefix, ext);
            if std::path::Path::new(&path).exists() {
                eprintln!("Loading background filter: {}", path);
                if let Ok(filter) = DomainBloomFilter::load(&path) {
                    bg_filters.push((filter, taxid));
                } else {
                    eprintln!("Failed to load {}", path);
                }
            }
        }
    }

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

    let scanner = MinimizerScanner::new(db.k(),db.l(),db.spaced_seed_mask,db.toggle_mask);
    let classify_start = std::time::Instant::now();
    let num_samples = config.samples.len();

    // Takes in either single FASTQ file or a dir of FASTQ files
    for sample in config.samples {
        let is_paired = sample.r2.is_some();
        let r1_path = sample.r1.to_string_lossy().to_string();
        let r2_path = sample.r2.as_ref().map(|p| p.to_string_lossy().to_string());   
        
        eprint!("\n>>> Processing sample ");
        let out_prefix = if num_samples == 1 && !std::path::Path::new(&config.output_prefix).is_dir() {
            eprint!("...");
            config.output_prefix.clone()
        } else {
            eprint!("{}", sample.name);
            format!("{}/{}", config.output_prefix, sample.name)
        };

        eprintln!();
        let no_output = config.no_output;
        let out_prefix_clone = out_prefix.clone();

        // thread only to write, no lock, so less overhead
        // we need sync channels here, since rayon worker threads will be much faster and more, just spitting
        // BatchResults into writer thread, but if we dont have sync channels, rayon threads will just stall...
        let (writer_tx, writer_rx) = std::sync::mpsc::sync_channel::<BatchResult>(config.threads * 4);
        // writer thread will act as a mutex, but without contention
        let writer_handle = std::thread::spawn(move || -> Result<(HashMap<u32, usize>, usize, usize, usize, usize)> {
            //now dont have to worry about flushing it out, since BufWriter will automatically write to disk and reset
            //when it;s full
            let mut writer = if !no_output {
                let output_file = File::create(format!("{}.output", out_prefix_clone))?;
                Some(io::BufWriter::with_capacity(4 * 1024 * 1024, output_file))
            } else {
                None
            };            
            
            let mut global_report: HashMap<u32, usize> = HashMap::new();
            // I could use AtomicUsize for these, since these are single values, but too lazy
            let mut total_c: usize = 0;
            let mut total_u: usize = 0;
            let mut total_a: usize = 0;
            let mut total_n: usize = 0;

            for batch in writer_rx {
                if let Some(ref mut w) = writer {
                    w.write_all(&batch.output_data)?;
                }
                for (taxid, count) in batch.report_counts {
                    *global_report.entry(taxid).or_default() += count;
                }
                total_c += batch.classified;
                total_u += batch.unclassified;
                total_a += batch.ambiguous;
                total_n += batch.num_records;
            }

            if let Some(ref mut w) = writer {
                w.flush()?;
            }
            Ok((global_report, total_c, total_u, total_a, total_n))
        });

        let batch_size = 5_000;
        // by doing this, reader keep pushes batches into crossbeam channel, and bridge lets each rayon worker
        // grab one batch whenever it's free (mpmc), so no contention like before, where if one thread slow, then all had 
        // to wait before starting on the nextc chunk (mpsc)
        let (batch_tx, batch_rx) = crossbeam_channel::bounded::<ReadBatch>(config.threads * 2);

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
                config.minimum_hit_groups,
                is_paired,
                &bg_filters,
                config.no_output
            );
            let _ = writer_tx.send(result);
        });

        drop(writer_tx);
        let (global_report, c, u, a, num_records) = match writer_handle.join() {
            Ok(Ok(v)) => v,
            Ok(Err(e)) => return Err(e),
            Err(_) => anyhow::bail!("Writer thread panicked"),
        };

        match reader_handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e),
            Err(_) => anyhow::bail!("Reader thread panicked"),
        }

        {
            let f = File::create(format!("{}.report", out_prefix))?;
            let mut writer = io::BufWriter::new(f);
            let mut sorted: Vec<(u32, usize)> = global_report.into_iter().collect();
            sorted.sort_by(|a, b| b.1.cmp(&a.1));

            writeln!(writer, "Target_TaxID\tName\tRead_Count\tRatio")?;
            for (taxid, count) in &sorted {
                let ratio = if num_records > 0 { *count as f64 / num_records as f64 } else { 0.0 };
                let taxid_label = if *taxid == u32::MAX { "Ambiguous".to_string() } else { taxid.to_string() };
                let name = db.get_taxid_name(*taxid);
                
                writeln!(writer, "{}\t{}\t{}\t{:0.3}", taxid_label, name, count, ratio)?;
            }
            writer.flush()?;
        }

        if c + u + a != num_records {
            eprintln!("WARNING: C+U+A is not adding up to number of reads");
        }
        eprintln!("Classified:   {}", c);
        eprintln!("Ambiguous:    {}", a);
        eprintln!("Unclassified: {}", u);
        eprintln!("Total:        {}", num_records);
        eprintln!();
    }
    let classify_elapsed = classify_start.elapsed();

    eprintln!("Took {:.2}s to classify sample(s)", classify_elapsed.as_secs_f64());
    Ok(())
}

#[inline(always)]
fn apply_penalty(hit_groups: usize, penalty: usize) -> usize {
    if hit_groups > penalty { hit_groups - penalty } else { 0 }
}

// Tried not to use mutex here since it was keep giving me extra lock contention overhead
// This was the main reason I switched to having a main writer thread
#[inline(never)]
fn classify_batch(
    batch: &ReadBatch,
    db: &KmerDatabase,
    scanner: &MinimizerScanner,
    acc_registry: &Option<AccessionRegistry>,
    min_hit_groups: usize,
    is_paired: bool,
    bg_filters: &[(DomainBloomFilter, u32)],
    no_output: bool,
) -> BatchResult {
    let batch_len = batch.batch1.offsets.len();    
    let mut local_output: Vec<u8> = if no_output { Vec::new() } else { Vec::with_capacity(batch_len * 200) };
    let mut local_report: FxHashMap<u32, usize> = FxHashMap::default();
    let mut local_classified: usize = 0;
    let mut local_unclassified: usize = 0;
    let mut local_ambiguous: usize = 0;

    let mut hits1: Vec<Hit> = Vec::new();
    let mut hits2: Vec<Hit> = Vec::new();
    let mut minimizer_buf1: Vec<u64> = Vec::new();
    let mut minimizer_buf2: Vec<u64> = Vec::new();

    for i in 0..batch_len {
        let (h1_start, h1_end, s1_start, s1_end) = batch.batch1.offsets[i];
        let header1 = &batch.batch1.data[h1_start..h1_end];
        let seq1 = &batch.batch1.data[s1_start..s1_end];

        hits1.clear();
        db.query_into(scanner, seq1, &mut minimizer_buf1, &mut hits1, bg_filters);

        let (seq2_len, has_r2) = if let Some(ref b2) = batch.batch2 {
            let (_, _, s2_start, s2_end) = b2.offsets[i];
            let seq2 = &b2.data[s2_start..s2_end];
            hits2.clear();
            db.query_into(scanner, seq2, &mut minimizer_buf2, &mut hits2, bg_filters);
            (seq2.len(), true)
        } else {
            (0, false)
        };

        let write_prefix = |out: &mut Vec<u8>, status: &str, taxid: &dyn std::fmt::Display| {
            let header_str = unsafe { std::str::from_utf8_unchecked(&header1) };
            if is_paired {
                let _ = write!(out, "{}\t{}\t{}\t{}|{}\t", status, header_str, taxid, seq1.len(), seq2_len);
            } else {
                let _ = write!(out, "{}\t{}\t{}\t{}\t", status, header_str, taxid, seq1.len());
            }
        };

        let mut valid_hits: usize = 0;
        let mut taxid_counts = [0u32; 256];
        let mut bg_counts = FxHashMap::default();
        let mut target_hit_groups: usize = 0;
        let mut bg_hit_groups: usize = 0;

        // Inline closure to process hits without code duplication
        let mut process_hits = |hits: &[Hit], minimizer_buf: &[u64]| {
            let mut last_minimizer: u64 = u64::MAX;
            for (idx, hit) in hits.iter().enumerate() {
                if hit.is_hit {
                    valid_hits += 1;
                    taxid_counts[hit.taxid_idx as usize] += 1;
                    if idx < minimizer_buf.len() {
                        let m = minimizer_buf[idx];
                        if m != last_minimizer {
                            target_hit_groups += 1;
                            last_minimizer = m;
                        }
                    }
                } else if hit.bg_taxid != 0 {
                    valid_hits += 1;
                    *bg_counts.entry(hit.bg_taxid).or_default() += 1;
                    if idx < minimizer_buf.len() {
                        let m = minimizer_buf[idx];
                        if m != last_minimizer {
                            bg_hit_groups += 1;
                            last_minimizer = m;
                        }
                    }
                }
            }
        };
        process_hits(&hits1, &minimizer_buf1); // this is for r1
        if has_r2 {
            process_hits(&hits2, &minimizer_buf2); // this is for r2
        }

        let effective_hit_groups = target_hit_groups + apply_penalty(bg_hit_groups, 1);

        //Unclassified
        if valid_hits == 0 || effective_hit_groups < min_hit_groups {
            if !no_output {
                write_prefix(&mut local_output, "U", &0);
                write_hit_pattern(&mut local_output, &hits1, db, &acc_registry);
                if has_r2 {
                    local_output.extend_from_slice(b" |:| ");
                    write_hit_pattern(&mut local_output, &hits2, db, &acc_registry);
                }
                local_output.push(b'\n');
            }
            local_unclassified += 1;
        } else {
            let n = db.index_to_taxid.len();
            
            let mut active_indices = Vec::new();
            for (idx, &count) in taxid_counts.iter().enumerate() {
                if count > 0 && idx < n {
                    active_indices.push(idx as u8);
                }
            }

            let mut path_weights = vec![0u32; n];
            let mut max_target_weight = 0;

            for i in 0..n {
                let mut w = 0;
                for &j in &active_indices {
                    // If j is an ancestor of i or vice versa, hit contribute to i's path
                    if db.is_ancestor(j, i as u8) {
                        w += taxid_counts[j as usize];
                    }
                }
                path_weights[i] = w;
                if w > max_target_weight {
                    max_target_weight = w;
                }
            }

            let mut best_target_nodes = Vec::new();
            if max_target_weight > 0 {
                for i in 0..n {
                    if path_weights[i] == max_target_weight {
                        best_target_nodes.push(i as u8);
                    }
                }
            }

            let mut best_bg_taxids: Vec<u32> = Vec::new();
            let mut best_bg_count: u32 = 0;
            for (&taxid, &count) in &bg_counts {
                if count > best_bg_count {
                    best_bg_count = count;
                    best_bg_taxids.clear();
                    best_bg_taxids.push(taxid);
                } else if count == best_bg_count {
                    best_bg_taxids.push(taxid);
                }
            }

            // tie breaker
            if max_target_weight >= best_bg_count && max_target_weight > 0 {
                
                let mut resolved_taxid_idx: Option<u8> = None;
                
                if best_target_nodes.len() == 1 {
                    resolved_taxid_idx = Some(best_target_nodes[0]);
                } else if best_target_nodes.len() > 1 {
                    // Doing LCA here, but LCA must be an ancestor all tied nodes
                    let mut candidate_lcas = Vec::new();
                    for k in 0..n {
                        let mut is_common = true;
                        for &s in &best_target_nodes {
                            if !db.is_ancestor(k as u8, s) {
                                is_common = false;
                                break;
                            }
                        }
                        if is_common {
                            candidate_lcas.push(k as u8);
                        }
                    }

                    if !candidate_lcas.is_empty() {
                        let mut lca_idx = candidate_lcas[0];
                        for &c in &candidate_lcas[1..] {
                            if db.is_ancestor(lca_idx, c) {
                                lca_idx = c; // c is deeper
                            }
                        }
                        resolved_taxid_idx = Some(lca_idx);
                    }
                }

                // Classified
                if let Some(final_idx) = resolved_taxid_idx {
                    let best_taxid = db.true_taxid(final_idx);
                    if !no_output {
                        write_prefix(&mut local_output, "C", &best_taxid);
                        write_hit_pattern(&mut local_output, &hits1, db, &acc_registry);
                        if has_r2 {
                            local_output.extend_from_slice(b" |:| ");
                            write_hit_pattern(&mut local_output, &hits2, db, &acc_registry);
                        }
                        local_output.push(b'\n');
                    }
                    *local_report.entry(best_taxid).or_default() += 1;
                    local_classified += 1;
                } else {
                    // no LCA then it's ambiguous
                    if !no_output {
                        let tied_taxids_str = best_target_nodes.iter().map(|&idx| db.true_taxid(idx).to_string()).join(",");
                        write_prefix(&mut local_output, "A", &tied_taxids_str);
                        write_hit_pattern(&mut local_output, &hits1, db, &acc_registry);
                        if has_r2 {
                            local_output.extend_from_slice(b" |:| ");
                            write_hit_pattern(&mut local_output, &hits2, db, &acc_registry);
                        }
                        local_output.push(b'\n');
                    }
                    *local_report.entry(u32::MAX).or_default() += 1; 
                    local_ambiguous += 1;
                }
                
            } else if best_bg_count > 0 {
                let final_bg_taxid = if best_bg_taxids.len() == 1 {
                    best_bg_taxids[0]
                } else {
                    resolve_domain_lca(best_bg_taxids.iter())
                };

                // background classified
                if !no_output {
                    write_prefix(&mut local_output, "C", &final_bg_taxid);
                    write_hit_pattern(&mut local_output, &hits1, db, &acc_registry);
                    if has_r2 {
                        local_output.extend_from_slice(b" |:| ");
                        write_hit_pattern(&mut local_output, &hits2, db, &acc_registry);
                    }
                    local_output.push(b'\n');
                }
                
                *local_report.entry(final_bg_taxid).or_default() += 1;
                local_classified += 1;
            }
        }
    }

    BatchResult {
        output_data: local_output,
        report_counts: local_report,
        classified: local_classified,
        unclassified: local_unclassified,
        ambiguous: local_ambiguous,
        num_records: batch_len,
    }
}

fn write_hit_pattern(out: &mut Vec<u8>, hits: &[Hit], db: &KmerDatabase, acc_registry: &Option<AccessionRegistry>) {
    if hits.is_empty() {
        out.extend_from_slice(b"0:0");
        return;
    }

    let get_taxid = |hit: &Hit| -> u32 {
        if hit.is_hit { db.true_taxid(hit.taxid_idx) }
        else if hit.bg_taxid != 0 { hit.bg_taxid }
        else { 0 }
    };

    let mut current_taxid = get_taxid(&hits[0]);
    let mut run_len: u32 = 1;

    let mut current_class_id = hits[0].accession_class_id;
    let reg = acc_registry.as_ref();

    for hit in hits.iter().skip(1) {
        let t = get_taxid(hit);
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

fn read_single_batches(r1_path: &str, batch_size: usize, tx: &crossbeam_channel::Sender<ReadBatch>) -> Result<()> {
    let mut r1_reader = open_fastx(r1_path)?;
    loop {
        let recs1 = read_flat_batch(&mut r1_reader, batch_size)?;
        if recs1.offsets.is_empty() { break; }
        if tx.send(ReadBatch { batch1: recs1, batch2: None }).is_err() { break; }
    }
    Ok(())
}

fn read_paired_batches(r1_path: &str, r2_path: &str, batch_size: usize, tx: &crossbeam_channel::Sender<ReadBatch>) -> Result<()> {
    let (r2_tx, r2_rx) = crossbeam_channel::bounded::<Option<FlatBatch>>(4);
    let r2_path_owned = r2_path.to_string();
    let bs = batch_size;

    let r2_handle = std::thread::spawn(move || -> Result<()> {
        let mut r2_reader = open_fastx(&r2_path_owned)?;
        loop {
            let recs2 = read_flat_batch(&mut r2_reader, bs)?;
            if recs2.offsets.is_empty() {
                let _ = r2_tx.send(None);
                break;
            }
            if r2_tx.send(Some(recs2)).is_err() { break; }
        }
        Ok(())
    });

    let mut r1_reader = open_fastx(r1_path)?;
    loop {
        let recs1 = read_flat_batch(&mut r1_reader, batch_size)?;
        if recs1.offsets.is_empty() { break; }
        let recs2 = match r2_rx.recv() {
            Ok(Some(r2)) => r2,
            Ok(None) => anyhow::bail!("R2 file ended before R1"),
            Err(_) => anyhow::bail!("R2 reader thread disconnected"),
        };

        if recs2.offsets.len() != recs1.offsets.len() {
            anyhow::bail!("There is a read count mismatch of paired end reads");
        }
        if tx.send(ReadBatch { batch1: recs1, batch2: Some(recs2) }).is_err() { break; }
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
fn read_flat_batch(reader: &mut FastxReader, n: usize) -> Result<FlatBatch> {
    
    let mut data = Vec::with_capacity(n * 250);
    let mut offsets = Vec::with_capacity(n);

    for _ in 0..n {
        match reader.next() {
            Some(Ok(rec)) => {
                let seq = rec.seq();
                if !seq.is_empty() {
                    let h_start = data.len();
                    data.extend_from_slice(rec.id());
                    let h_end = data.len();

                    let s_start = data.len();
                    data.extend_from_slice(&seq);
                    let s_end = data.len();

                    offsets.push((h_start, h_end, s_start, s_end));
                }
            }
            Some(Err(e)) => return Err(anyhow::anyhow!("Error reading records: {}", e)),
            None => break,
        }
    }
    
    Ok(FlatBatch { data, offsets })
}