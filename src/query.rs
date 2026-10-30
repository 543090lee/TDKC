pub mod read_finder;
use crate::database::TaxIdx;
use std::io::{self, Write, BufReader};
use std::fs::File;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use anyhow::Result;
use rayon::prelude::*;
use crate::database::{AccessionRegistry, DomainBloomFilter};
use crate::querydb::{Lookup, QueryDb};
use crate::minimizer::MinimizerScanner;
use rustc_hash::FxHashMap;
use crate::utils::init_thread_pool;
use crate::taxonomy::{TAXID_ARCHAEA, TAXID_BACTERIA, TAXID_FUNGI, TAXID_VIRAL};
use read_finder::Sample;
use std::process::{Command, Stdio};

pub struct QueryConfig {
    pub db_dir: String,
    pub samples: Vec<Sample>,
    pub threads: usize,
    pub use_accessions: bool,
    pub background: bool,
    pub minimum_hit_groups: usize,
    pub output_prefix: String,
    pub no_output: bool,
    pub gzip_output: bool,
    pub gzip_multithreaded: bool,
    pub domain_penalty: usize,
    pub counts: bool,
}

pub struct FlatBatch {
    pub data: Vec<u8>,
    pub offsets: Vec<(usize, usize, usize, usize)>,
}

pub struct ReadBatch {
    pub batch1: FlatBatch,
    pub batch2: Option<FlatBatch>,
}

type AccCountMap = FxHashMap<u32, (u32, u32)>;
struct SeenBitset {
    words: Vec<AtomicU64>,
}

impl SeenBitset {
    fn new(num_dense: usize) -> Self {
        let num_words = (num_dense + 63) / 64;
        let mut words = Vec::with_capacity(num_words);
        for _ in 0..num_words {
            words.push(AtomicU64::new(0));
        }
        Self { words }
    }

    // Set the bit then returns true iff this thread was the one that flipped 0 to 1
    #[inline]
    fn test_and_set(&self, dense_pos: u32) -> bool {
        let word = (dense_pos as usize) / 64;
        let bit = 1u64 << ((dense_pos as u64) % 64);
        let prev = self.words[word].fetch_or(bit, Ordering::Relaxed);
        prev & bit == 0
    }
}

// this goes to writer thread
struct BatchResult {
    output_data: Vec<u8>,
    report_counts: Vec<usize>,
    acc_counts: AccCountMap,
    classified: usize,
    unclassified: usize,
    ambiguous: usize,
    num_records: usize,
}

const GROUP: usize = 8;

pub fn run_query(config: QueryConfig) -> Result<()> {
    init_thread_pool(config.threads);

    let db_prefix = format!("{}/db", config.db_dir);
    let db = QueryDb::load(&db_prefix, config.use_accessions)?;

    let mut bg_filters = Vec::new();
    if config.background {
        let domains = [
            ("bacteria.bloom", TAXID_BACTERIA),
            ("archaea.bloom", TAXID_ARCHAEA),
            ("viral.bloom", TAXID_VIRAL),
            ("fungi.bloom", TAXID_FUNGI)
        ];
        for (ext, taxid) in domains {
            let path = format!("{}.{}", db_prefix, ext);
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

    let acc_path = format!("{}.accessions", db_prefix);
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

    if config.counts && (acc_registry.is_none() || db.accessions.is_none()) {
        anyhow::bail!("--counts requires an accession database");
    }

    let scanner = MinimizerScanner::new(db.k, db.l, db.spaced_seed_mask, db.toggle_mask);
    let classify_start = std::time::Instant::now();
    let num_samples = config.samples.len();
    let n_taxa = db.index_to_taxid.len();

    let ctx = ClassifyCtx {
        db: &db,
        scanner: &scanner,
        acc_registry: acc_registry.as_ref(),
        min_hit_groups: config.minimum_hit_groups,
        bg_filters: &bg_filters,
        no_output: config.no_output,
        gzip_output: config.gzip_output,
        domain_penalty: config.domain_penalty,
    };

    // Takes in either single FASTQ file or a dir of FASTQ files
    for sample in &config.samples {
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
        let output_ext = if config.gzip_output { ".output.gz" } else { ".output" };
        let out_prefix_clone = out_prefix.clone();

        let seen_bits: Option<Arc<SeenBitset>> = if config.counts {
            db.accessions.as_ref().map(|acc| Arc::new(SeenBitset::new(acc.num_dense())))
        } else {
            None
        };

        // thread only to write, no lock, so less overhead
        // we need sync channels here, since rayon worker threads will be much faster and more, just spitting
        // BatchResults into writer thread, but if we dont have sync channels, rayon threads will just stall...
        let (writer_tx, writer_rx) = std::sync::mpsc::sync_channel::<BatchResult>(config.threads * 4);
        let want_counts = config.counts;
        // writer thread will act as a mutex, but without contention
        let writer_handle = std::thread::spawn(move || -> Result<(Vec<usize>, AccCountMap, usize, usize, usize, usize)> {
            //now dont have to worry about flushing it out, since BufWriter will automatically write to disk and reset
            //when it;s full
            let mut writer = if !no_output {
                let output_file = File::create(format!("{}{}", out_prefix_clone, output_ext))?;
                Some(io::BufWriter::with_capacity(4 * 1024 * 1024, output_file))
            } else {
                None
            };

            let mut global_report: Vec<usize> = vec![0; n_taxa + 1];
            let mut global_acc_counts: AccCountMap = FxHashMap::default();
            // I could use AtomicUsize for these, since these are single values, but too lazy
            let mut total_c: usize = 0;
            let mut total_u: usize = 0;
            let mut total_a: usize = 0;
            let mut total_n: usize = 0;

            for batch in writer_rx {
                if let Some(ref mut w) = writer {
                    w.write_all(&batch.output_data)?;
                }
                for (slot, count) in global_report.iter_mut().zip(batch.report_counts.iter()) {
                    *slot += *count;
                }
                if want_counts {
                    for (acc_id, (total, uniq)) in batch.acc_counts {
                        let slot = global_acc_counts.entry(acc_id).or_insert((0, 0));
                        slot.0 += total;
                        slot.1 += uniq;
                    }
                }
                total_c += batch.classified;
                total_u += batch.unclassified;
                total_a += batch.ambiguous;
                total_n += batch.num_records;
            }

            if let Some(ref mut w) = writer {
                w.flush()?;
            }
            Ok((global_report, global_acc_counts, total_c, total_u, total_a, total_n))
        });

        let batch_size = 5_000;
        // by doing this, reader keep pushes batches into crossbeam channel, and bridge lets each rayon worker
        // grab one batch whenever it's free (mpmc), so no contention like before, where if one thread slow, then all had
        // to wait before starting on the nextc chunk (mpsc)
        let (batch_tx, batch_rx) = crossbeam_channel::bounded::<ReadBatch>(config.threads * 2);

        let gz_threads = if is_paired { (config.threads / 2).clamp(1, 4) } else { config.threads.clamp(1, 4) };
        let gzip_multithreaded = config.gzip_multithreaded;
        let reader_handle = std::thread::spawn(move || -> Result<()> {
            if let Some(ref r2_path_str) = r2_path {
                read_paired_batches(&r1_path, r2_path_str, batch_size, &batch_tx, gzip_multithreaded, gz_threads)?;
            } else {
                read_single_batches(&r1_path, batch_size, &batch_tx, gzip_multithreaded, gz_threads)?;
            }
            // sender is done, so done reading (EOF)
            drop(batch_tx);
            Ok(())
        });

        //par_bridge() let's channel to give batch to rayon thread whenever it's free
        let seen_ref = seen_bits.as_deref();
        batch_rx.into_iter().par_bridge().for_each_init(
            || Scratch::new(n_taxa),
            |scratch, batch| {
                let result = classify_batch(&batch, &ctx, is_paired, seen_ref, scratch);
                let _ = writer_tx.send(result);
            },
        );

        drop(writer_tx);
        let (global_report, global_acc_counts, c, u, a, num_records) = match writer_handle.join() {
            Ok(Ok(v)) => v,
            Ok(Err(e)) => return Err(e),
            Err(_) => anyhow::bail!("Writer thread panicked"),
        };

        match reader_handle.join() {
            Ok(Ok(())) => {}
            Ok(Err(e)) => return Err(e),
            Err(_) => anyhow::bail!("Reader thread panicked"),
        }

        // .report
        {
            let f = File::create(format!("{}.report", out_prefix))?;
            let mut writer = io::BufWriter::new(f);
            let mut sorted: Vec<(u32, usize)> = global_report
                .iter()
                .enumerate()
                .filter(|(_, &c)| c > 0)
                .map(|(i, &c)| (if i == n_taxa { u32::MAX } else { db.index_to_taxid[i] }, c))
                .collect();
            sorted.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

            writeln!(writer, "Target_TaxID\tName\tRead_Count\tRatio")?;
            for (taxid, count) in &sorted {
                let ratio = if num_records > 0 { *count as f64 / num_records as f64 } else { 0.0 };
                let taxid_label = if *taxid == u32::MAX { "Ambiguous".to_string() } else { taxid.to_string() };
                let name = db.get_taxid_name(*taxid);

                writeln!(writer, "{}\t{}\t{}\t{:0.3}", taxid_label, name, count, ratio)?;
            }
            writer.flush()?;
        }

        // strain.txt
        if config.counts {
            if let Some(reg) = acc_registry.as_ref() {
                let strain_path = format!("{}.strain.txt", out_prefix);
                let f = File::create(&strain_path)?;
                let mut writer = io::BufWriter::new(f);
                writeln!(writer, "Accession\tTotal_Hits\tUnique_Hits")?;

                let mut rows: Vec<(u32, u32, u32)> = global_acc_counts
                    .into_iter()
                    .map(|(acc_id, (total, uniq))| (acc_id, total, uniq))
                    .collect();
                rows.sort_by(|a, b| {
                    b.1.cmp(&a.1).then(b.2.cmp(&a.2)).then(a.0.cmp(&b.0))
                });

                for (acc_id, total, uniq) in &rows {
                    writeln!(writer, "{}\t{}\t{}", reg.get_name(*acc_id), total, uniq)?;
                }
                writer.flush()?;
            }
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

struct ClassifyCtx<'a> {
    db: &'a QueryDb,
    scanner: &'a MinimizerScanner,
    acc_registry: Option<&'a AccessionRegistry>,
    min_hit_groups: usize,
    bg_filters: &'a [(DomainBloomFilter, u32)],
    no_output: bool,
    gzip_output: bool,
    domain_penalty: usize,
}

struct Scratch {
    minimizer_buf: Vec<u64>,
    runs: Vec<(u64, u32)>,
    run_lookup: Vec<u32>,
    lookups: Vec<Lookup>,
    lk_scratch: Vec<u32>,
    spans: Vec<(u32, u32, u32, u32)>,
    taxid_counts: Vec<u32>,
    path_weights: Vec<u32>,
    best: Vec<TaxIdx>,
    cand: Vec<u64>,
    acc_buf: Vec<u32>,
    itoa: itoa::Buffer,
}

impl Scratch {
    fn new(n_taxa: usize) -> Self {
        Self {
            minimizer_buf: Vec::with_capacity(512),
            runs: Vec::with_capacity(GROUP * 128),
            run_lookup: Vec::with_capacity(GROUP * 128),
            lookups: Vec::with_capacity(GROUP * 128),
            lk_scratch: Vec::with_capacity(GROUP * 128),
            spans: Vec::with_capacity(GROUP),
            taxid_counts: vec![0u32; n_taxa.max(1)],
            path_weights: vec![0u32; n_taxa.max(1)],
            best: Vec::with_capacity(n_taxa.max(1)),
            cand: vec![0u64; (n_taxa.max(1) + 63) / 64],
            acc_buf: Vec::new(),
            itoa: itoa::Buffer::new(),
        }
    }
}

#[inline(always)]
fn scan_append(scanner: &MinimizerScanner, seq: &[u8], sc: &mut Scratch) {
    let mut buf = std::mem::take(&mut sc.minimizer_buf);
    scanner.scan_into(seq, &mut buf);
    let mut cur: u64 = 0;
    let mut len: u32 = 0;
    for &m in &buf {
        if len > 0 && m == cur {
            len += 1;
            continue;
        }
        if len > 0 {
            push_run(sc, cur, len);
        }
        cur = m;
        len = 1;
    }
    if len > 0 {
        push_run(sc, cur, len);
    }
    sc.minimizer_buf = buf;
}

#[inline(always)]
fn push_run(sc: &mut Scratch, m: u64, len: u32) {
    sc.runs.push((m, len));
    if m == u64::MAX {
        sc.run_lookup.push(u32::MAX);
    } else {
        sc.run_lookup.push(sc.lookups.len() as u32);
        sc.lookups.push(Lookup::new(m));
    }
}

// Tried not to use mutex here since it was keep giving me extra lock contention overhead
// This was the main reason I switched to having a main writer thread
#[inline(never)]
fn classify_batch(
    batch: &ReadBatch,
    ctx: &ClassifyCtx,
    is_paired: bool,
    seen_bits: Option<&SeenBitset>,
    sc: &mut Scratch,
) -> BatchResult {
    let db = ctx.db;
    let batch_len = batch.batch1.offsets.len();
    let n = db.index_to_taxid.len();
    let mut local_output: Vec<u8> = if ctx.no_output { Vec::new() } else { Vec::with_capacity(batch_len * 160) };
    let mut local_report: Vec<usize> = vec![0; n + 1];
    let mut local_acc_counts: AccCountMap = FxHashMap::default();
    let mut local_classified: usize = 0;
    let mut local_unclassified: usize = 0;
    let mut local_ambiguous: usize = 0;

    let acc_for_counts = if seen_bits.is_some() { db.accessions.as_ref() } else { None };
    let b1 = &batch.batch1;
    let b2 = batch.batch2.as_ref();

    let mut g0 = 0usize;
    while g0 < batch_len {
        let g1 = (g0 + GROUP).min(batch_len);

        sc.runs.clear();
        sc.run_lookup.clear();
        sc.lookups.clear();
        sc.spans.clear();
        for i in g0..g1 {
            let (_, _, s1_start, s1_end) = b1.offsets[i];
            let r1s = sc.runs.len() as u32;
            scan_append(ctx.scanner, &b1.data[s1_start..s1_end], sc);
            let r1e = sc.runs.len() as u32;
            let (r2s, r2e) = if let Some(b2) = b2 {
                let (_, _, s2_start, s2_end) = b2.offsets[i];
                let a = sc.runs.len() as u32;
                scan_append(ctx.scanner, &b2.data[s2_start..s2_end], sc);
                (a, sc.runs.len() as u32)
            } else {
                (r1e, r1e)
            };
            sc.spans.push((r1s, r1e, r2s, r2e));
        }
        db.lookup_batch(&mut sc.lookups, &mut sc.lk_scratch, ctx.bg_filters);

        for (gi, i) in (g0..g1).enumerate() {
            let (r1s, r1e, r2s, r2e) = sc.spans[gi];
            let (h1_start, h1_end, s1_start, s1_end) = b1.offsets[i];
            let header = &b1.data[h1_start..h1_end];
            let seq1_len = s1_end - s1_start;
            let seq2_len = b2.map_or(0, |b| { let (_, _, s, e) = b.offsets[i]; e - s });

            if let (Some(sb), Some(acc)) = (seen_bits, acc_for_counts) {
                for r in (r1s as usize..r1e as usize).chain(r2s as usize..r2e as usize) {
                    let li = sc.run_lookup[r];
                    if li == u32::MAX {
                        continue;
                    }
                    let lk = &sc.lookups[li as usize];
                    // Helper that credits an accession-class for the very first time seen
                    if lk.is_hit && lk.class_id != u32::MAX && sb.test_and_set(lk.dense_pos) {
                        acc.decode_class_into(lk.class_id, &mut sc.acc_buf);
                        let is_unique = sc.acc_buf.len() == 1;
                        for &a in &sc.acc_buf {
                            let slot = local_acc_counts.entry(a).or_insert((0, 0));
                            slot.0 = slot.0.saturating_add(1);
                            if is_unique {
                                slot.1 = slot.1.saturating_add(1);
                            }
                        }
                    }
                }
            }

            let mut valid_hits: usize = 0;
            let counts = &mut sc.taxid_counts[..];
            counts[..n].fill(0);
            let mut target_hit_groups: usize = 0;
            let mut bg_hit_groups: usize = 0;

            for (rs, re) in [(r1s, r1e), (r2s, r2e)] {
                let mut last_minimizer: u64 = u64::MAX;
                for r in rs as usize..re as usize {
                    let li = sc.run_lookup[r];
                    if li == u32::MAX {
                        continue;
                    }
                    let (m, len) = sc.runs[r];
                    let lk = &sc.lookups[li as usize];
                    if lk.is_hit {
                        valid_hits += len as usize;
                        let ti = lk.taxid_idx as usize;
                        if ti < n {
                            counts[ti] += len;
                        }
                        if m != last_minimizer {
                            target_hit_groups += 1;
                            last_minimizer = m;
                        }
                    } else if lk.bg_taxid != 0 {
                        valid_hits += len as usize;
                        if let Some(di) = db.get_taxid_index(lk.bg_taxid) {
                            counts[di as usize] += len;
                        }
                        if m != last_minimizer {
                            bg_hit_groups += 1;
                            last_minimizer = m;
                        }
                    }
                }
            }

            let effective_hit_groups = target_hit_groups + apply_penalty(bg_hit_groups, ctx.domain_penalty);

            //Unclassified
            if valid_hits == 0 || effective_hit_groups < ctx.min_hit_groups {
                if !ctx.no_output {
                    write_prefix(&mut local_output, b"U", header);
                    local_output.extend_from_slice(b"0");
                    write_patterns(&mut local_output, sc, r1s, r1e, r2s, r2e, b2.is_some(), ctx, is_paired, seq1_len, seq2_len);
                }
                local_unclassified += 1;
                continue;
            }

            let pw = &mut sc.path_weights[..];
            pw[..n].fill(0);
            let mut max_w = 0u32;
            for j in 0..n {
                let cj = counts[j];
                if cj > 0 {
                    // If j is an ancestor of i or vice versa, hit contribute to i's path
                    for &d in &db.desc_lists[j] {
                        let v = pw[d as usize] + cj;
                        pw[d as usize] = v;
                        if v > max_w {
                            max_w = v;
                        }
                    }
                }
            }

            if max_w == 0 {
                continue;
            }

            // tie breaker
            let best = &mut sc.best;
            best.clear();
            for i2 in 0..n {
                if pw[i2] == max_w {
                    best.push(i2 as TaxIdx);
                }
            }
            let nbest = best.len();

            let resolved: Option<TaxIdx> = if nbest == 1 {
                Some(best[0])
            } else {
                // Doing LCA here, but LCA must be an ancestor all tied nodes
                let nw = db.mask_words;
                let cand = &mut sc.cand[..nw];
                cand.fill(u64::MAX);
                for &s in best.iter() {
                    let m = db.anc_mask(s as usize);
                    for w in 0..nw {
                        cand[w] &= m[w];
                    }
                }
                if ctx.bg_filters.is_empty() {
                    for w in 0..nw {
                        cand[w] &= !db.domain_mask[w];
                    }
                }
                let mut lca: Option<TaxIdx> = None;
                for w in 0..nw {
                    let mut bits = cand[w];
                    while bits != 0 {
                        let c = (w * 64 + bits.trailing_zeros() as usize) as TaxIdx;
                        bits &= bits - 1;
                        match lca {
                            None => lca = Some(c),
                            Some(cur) => {
                                if db.is_ancestor(cur, c) {
                                    lca = Some(c);
                                }
                            }
                        }
                    }
                }
                lca
            };

            match resolved {
                // Classified
                Some(final_idx) => {
                    let best_taxid = db.true_taxid(final_idx);
                    if !ctx.no_output {
                        write_prefix(&mut local_output, b"C", header);
                        local_output.extend_from_slice(sc.itoa.format(best_taxid).as_bytes());
                        write_patterns(&mut local_output, sc, r1s, r1e, r2s, r2e, b2.is_some(), ctx, is_paired, seq1_len, seq2_len);
                    }
                    local_report[final_idx as usize] += 1;
                    local_classified += 1;
                }
                // no LCA then it's ambiguous
                None => {
                    if !ctx.no_output {
                        write_prefix(&mut local_output, b"A", header);
                        for (bi, &b) in best[..nbest].iter().enumerate() {
                            if bi > 0 {
                                local_output.push(b',');
                            }
                            local_output.extend_from_slice(sc.itoa.format(db.true_taxid(b)).as_bytes());
                        }
                        write_patterns(&mut local_output, sc, r1s, r1e, r2s, r2e, b2.is_some(), ctx, is_paired, seq1_len, seq2_len);
                    }
                    local_report[n] += 1;
                    local_ambiguous += 1;
                }
            }
        }
        g0 = g1;
    }

    if ctx.gzip_output && !local_output.is_empty() {
        let mut enc = flate2::write::GzEncoder::new(Vec::with_capacity(local_output.len() / 4 + 64), flate2::Compression::fast());
        std::io::Write::write_all(&mut enc, &local_output).expect("gzip of in-memory buffer cannot fail");
        local_output = enc.finish().expect("gzip of in-memory buffer cannot fail");
    }

    BatchResult {
        output_data: local_output,
        report_counts: local_report,
        acc_counts: local_acc_counts,
        classified: local_classified,
        unclassified: local_unclassified,
        ambiguous: local_ambiguous,
        num_records: batch_len,
    }
}

#[inline(always)]
fn write_prefix(out: &mut Vec<u8>, status: &[u8], header: &[u8]) {
    out.extend_from_slice(status);
    out.push(b'\t');
    out.extend_from_slice(header);
    out.push(b'\t');
}

#[inline(never)]
fn write_patterns(
    out: &mut Vec<u8>,
    sc: &mut Scratch,
    r1s: u32,
    r1e: u32,
    r2s: u32,
    r2e: u32,
    has_r2: bool,
    ctx: &ClassifyCtx,
    is_paired: bool,
    len1: usize,
    len2: usize,
) {
    out.push(b'\t');
    out.extend_from_slice(sc.itoa.format(len1).as_bytes());
    if is_paired {
        out.push(b'|');
        out.extend_from_slice(sc.itoa.format(len2).as_bytes());
    }
    out.push(b'\t');
    write_hit_pattern(out, sc, r1s, r1e, ctx);
    if has_r2 {
        out.extend_from_slice(b" |:| ");
        write_hit_pattern(out, sc, r2s, r2e, ctx);
    }
    out.push(b'\n');
}

#[inline(always)]
fn write_hit_pattern(out: &mut Vec<u8>, sc: &mut Scratch, rs: u32, re: u32, ctx: &ClassifyCtx) {
    if rs == re {
        out.extend_from_slice(b"0:0");
        return;
    }
    let db = ctx.db;
    let mut cur_taxid: u32 = 0;
    let mut cur_class: u32 = u32::MAX;
    let mut run_len: u32 = 0;
    let mut first = true;
    for r in rs as usize..re as usize {
        let (_, len) = sc.runs[r];
        let li = sc.run_lookup[r];
        let (t, cls) = if li == u32::MAX {
            (0u32, u32::MAX)
        } else {
            let lk = &sc.lookups[li as usize];
            if lk.is_hit {
                (db.true_taxid(lk.taxid_idx), lk.class_id)
            } else if lk.bg_taxid != 0 {
                (lk.bg_taxid, lk.class_id)
            } else {
                (0, lk.class_id)
            }
        };
        if first {
            cur_taxid = t;
            cur_class = cls;
            run_len = len;
            first = false;
        } else if t == cur_taxid {
            run_len += len;
        } else {
            emit_run(out, cur_taxid, run_len, cur_class, ctx, &mut sc.acc_buf, &mut sc.itoa);
            cur_taxid = t;
            cur_class = cls;
            run_len = len;
        }
    }
    emit_run(out, cur_taxid, run_len, cur_class, ctx, &mut sc.acc_buf, &mut sc.itoa);
}

#[inline(always)]
fn emit_run(out: &mut Vec<u8>, taxid: u32, run_len: u32, class_id: u32, ctx: &ClassifyCtx, acc_buf: &mut Vec<u32>, itoa: &mut itoa::Buffer) {
    if let Some(reg) = ctx.acc_registry {
        if taxid != 0 && class_id != u32::MAX {
            if let Some(acc) = ctx.db.accessions.as_ref() {
                acc.decode_class_into(class_id, acc_buf);
                if !acc_buf.is_empty() {
                    out.extend_from_slice(itoa.format(taxid).as_bytes());
                    out.extend_from_slice(b"-{");
                    for (i, id) in acc_buf.iter().enumerate() {
                        if i > 0 {
                            out.push(b',');
                        }
                        out.extend_from_slice(reg.get_name(*id).as_bytes());
                    }
                    out.extend_from_slice(b"}:");
                    out.extend_from_slice(itoa.format(run_len).as_bytes());
                    out.push(b' ');
                    return;
                }
            }
        }
    }
    out.extend_from_slice(itoa.format(taxid).as_bytes());
    out.push(b':');
    out.extend_from_slice(itoa.format(run_len).as_bytes());
    out.push(b' ');
}

fn read_single_batches(r1_path: &str, batch_size: usize, tx: &crossbeam_channel::Sender<ReadBatch>, gzip_multithreaded: bool, gz_threads: usize) -> Result<()> {
    let mut r1_reader = if !gzip_multithreaded {
        open_fastx(r1_path, gz_threads)?
    } else {
        open_fastx_fast(r1_path)?
    };

    loop {
        let recs1 = read_flat_batch(&mut r1_reader, batch_size)?;
        if recs1.offsets.is_empty() { break; }
        if tx.send(ReadBatch { batch1: recs1, batch2: None }).is_err() { break; }
    }
    Ok(())
}

fn read_paired_batches(r1_path: &str, r2_path: &str, batch_size: usize, tx: &crossbeam_channel::Sender<ReadBatch>, gzip_multithreaded: bool, gz_threads: usize) -> Result<()> {
    let (r2_tx, r2_rx) = crossbeam_channel::bounded::<Option<FlatBatch>>(4);
    let r2_path_owned = r2_path.to_string();
    let bs = batch_size;

    let r2_handle = std::thread::spawn(move || -> Result<()> {
        let mut r2_reader = if !gzip_multithreaded {
            open_fastx(&r2_path_owned, gz_threads)?
        } else {
            open_fastx_fast(&r2_path_owned)?
        };

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

    let mut r1_reader = if !gzip_multithreaded {
        open_fastx(r1_path, gz_threads)?
    } else {
        open_fastx_fast(r1_path)?
    };
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

fn open_fastx(path: &str, gz_threads: usize) -> Result<FastxReader> {
    if path.ends_with(".gz") {
        if let Some(pgz) = crate::pgz::ParallelGzReader::open(path, gz_threads)? {
            let buf_reader = BufReader::with_capacity(4 * 1024 * 1024, pgz);
            let reader = needletail::parse_fastx_reader(buf_reader)
                .map_err(|e| anyhow::anyhow!("Cannot parse reads file {}: {}", path, e))?;
            return Ok(reader);
        }
    }
    let reader = needletail::parse_fastx_file(path)
        .map_err(|e| anyhow::anyhow!("Cannot open reads file {}", e))?;
    Ok(reader)
}

fn open_fastx_fast(path: &str) -> Result<FastxReader> {
    if path.ends_with(".gz") {
        let mut child = Command::new("rapidgzip")
            .args(["-d", "-c", "-P", "3", path])
            .stdout(Stdio::piped())
            .spawn()
            .map_err(|e| anyhow::anyhow!("is rapidgzip installed?: {}", e))?;

        let stdout = child.stdout.take().unwrap();

        let buf_reader = BufReader::with_capacity(4 * 1024 * 1024, stdout);

        let reader = needletail::parse_fastx_reader(buf_reader)
            .map_err(|e| anyhow::anyhow!("Cannot parse the stream: {}", e))?;

        Ok(reader)
    } else {
        let reader = needletail::parse_fastx_file(path)
            .map_err(|e| anyhow::anyhow!("Cannot open reads file {}", e))?;
        Ok(reader)
    }
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
