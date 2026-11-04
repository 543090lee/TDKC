use std::fs::File;
use std::io::{BufWriter, Read};
use std::path::Path;
use std::sync::atomic::{AtomicU64, AtomicU8, Ordering};
use std::sync::Arc;
use std::time::Instant;
use anyhow::{Context, Result};
use hyperloglogplus::{HyperLogLog, HyperLogLogPlus};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use crate::database::DomainBloomFilter;
use crate::minimizer::MinimizerScanner;
use crate::taxonomy::{TaxonomyTree, TAXID_ARCHAEA, TAXID_BACTERIA, TAXID_FUNGI, TAXID_VIRAL};
use crate::utils::{extract_accession, init_thread_pool};
use crate::hash::{BuildFmix64Hasher, bloom_hash_pair};

const RAW_FPR: f64 = 0.0001;

const BACTERIA: usize = 0;
const ARCHAEA: usize = 1;
const VIRAL: usize = 2;
const FUNGI: usize = 3;
const OTHER: usize = 4;
const NGROUPS: usize = 5;
const GROUP_NAMES: [&str; NGROUPS] = ["bacteria", "archaea", "viral", "fungi", "other"];

struct AtomicBloom {
    data: Vec<AtomicU8>,
    num_hashes: u32,
    num_bits: u64,
}

impl AtomicBloom {
    fn new(expected_items: usize, fpr: f64) -> Self {
        let n = (expected_items.max(1)) as f64;
        let ln2 = std::f64::consts::LN_2;
        let m = (-n * fpr.ln() / (ln2 * ln2)).ceil() as u64;
        let k = ((m as f64 / n) * ln2).round() as u32;
        let k = k.max(1);

        let num_bytes = ((m + 7) / 8) as usize;
        let mut data = Vec::with_capacity(num_bytes);
        for _ in 0..num_bytes {
            data.push(AtomicU8::new(0));
        }

        eprintln!("  {} bits ({:.2} GiB), {} hashes for {} items @ {}% FPR", m, num_bytes as f64 / 1_073_741_824.0, k, expected_items, fpr * 100.0);

        Self { data, num_hashes: k, num_bits: m }
    }

    #[inline]
    fn insert(&self, item: u64) {
        let (h1, h2) = bloom_hash_pair(item);
        for i in 0..self.num_hashes {
            let bit_pos = (h1.wrapping_add((i as u64).wrapping_mul(h2))) % self.num_bits;
            let byte_idx = (bit_pos / 8) as usize;
            let bit_mask = 1u8 << (bit_pos % 8);
            self.data[byte_idx].fetch_or(bit_mask, Ordering::Relaxed);
        }
    }

    #[inline]
    fn contains(&self, item: u64) -> bool {
        let (h1, h2) = bloom_hash_pair(item);
        for i in 0..self.num_hashes {
            let bit_pos = (h1.wrapping_add((i as u64).wrapping_mul(h2))) % self.num_bits;
            let byte_idx = (bit_pos / 8) as usize;
            let bit_mask = 1u8 << (bit_pos % 8);
            if self.data[byte_idx].load(Ordering::Relaxed) & bit_mask == 0 {
                return false;
            }
        }
        true
    }

    fn into_filter(self) -> DomainBloomFilter {
        let (k, m) = (self.num_hashes, self.num_bits);
        DomainBloomFilter::new(self.data.into_iter().map(|a| a.into_inner()).collect(), k, m)
    }
}

pub struct BuildDomainConfig {
    pub db_dir: String,
    pub threads: usize,
    pub fpr: f64,
    pub bacteria: bool,
    pub archaea: bool,
    pub viral: bool,
    pub fungi: bool,
}

struct RefFile {
    path: String,
    fixed: Option<usize>,
}

type Hll = HyperLogLogPlus<u64, BuildFmix64Hasher>;

fn new_hll() -> Hll {
    HyperLogLogPlus::<u64, BuildFmix64Hasher>::new(16, BuildFmix64Hasher::default()).unwrap()
}

pub fn run_build_domain(config: BuildDomainConfig) -> Result<()> {
    init_thread_pool(config.threads);

    let total_start = Instant::now();
    let db_path = std::path::PathBuf::from(&config.db_dir);
    let db_prefix = db_path.join("db").to_string_lossy().into_owned();
    let genome_dir = db_path.join("genome");

    let (k, l, spaced_seed_mask, toggle_mask) = load_db_meta(&db_prefix)?;
    let scanner = MinimizerScanner::new(k, l, spaced_seed_mask, toggle_mask);

    let none_specified = !config.bacteria && !config.archaea && !config.viral && !config.fungi;
    let want = |flag: bool| flag || none_specified;
    let mut requested: Vec<usize> = Vec::new();
    for (g, flag) in [(BACTERIA, config.bacteria), (ARCHAEA, config.archaea), (VIRAL, config.viral), (FUNGI, config.fungi)] {
        if !want(flag) {
            continue;
        }
        if genome_dir.join(format!("{}.fna", GROUP_NAMES[g])).exists() {
            requested.push(g);
        } else {
            eprintln!("Can't find {}.fna", GROUP_NAMES[g]);
        }
    }
    if requested.is_empty() {
        eprintln!("No domains specified.");
        return Ok(());
    }

    let mut ref_files: Vec<RefFile> = Vec::new();
    let mut names: Vec<String> = std::fs::read_dir(&genome_dir)
        .with_context(|| format!("Can't read {}", genome_dir.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".fna") && !n.starts_with('_'))
        .collect();
    names.sort();
    for name in names {
        let stem = name.trim_end_matches(".fna");
        let fixed = (0..OTHER).find(|&g| GROUP_NAMES[g] == stem);
        ref_files.push(RefFile { path: genome_dir.join(&name).to_string_lossy().into_owned(), fixed });
    }
    eprintln!("Reference files:");
    for f in &ref_files {
        eprintln!("  {}{}", f.path, f.fixed.map(|g| format!("  (all records = {})", GROUP_NAMES[g])).unwrap_or_default());
    }

    let t = Instant::now();
    let taxonomy = TaxonomyTree::load(&db_path.join("taxonomy/nodes.dmp").to_string_lossy())?;
    let prelim = crate::prep::load_prelim_map(&db_path.join("prelim_map.txt").to_string_lossy())?;
    let group_of = Arc::new(GroupResolver::new(taxonomy, prelim));
    eprintln!("  [taxonomy] {:.1}s", t.elapsed().as_secs_f64());

    eprintln!("\nPass A: estimating unique minimizers per domain group over the reference...");
    let t = Instant::now();
    let rec_counts: Vec<AtomicU64> = (0..NGROUPS).map(|_| AtomicU64::new(0)).collect();
    let mut hlls = stream_reference(&ref_files, &group_of, &rec_counts, |batch| {
        let mut local: Vec<Hll> = (0..NGROUPS).map(|_| new_hll()).collect();
        let mut mins = Vec::new();
        for (g, seq) in batch {
            scanner.scan_into(seq, &mut mins);
            for &m in &mins {
                if m != u64::MAX {
                    local[*g].insert(&m);
                }
            }
        }
        local
    })?;
    let mut raw_counts = [0usize; NGROUPS];
    for g in 0..NGROUPS {
        raw_counts[g] = hlls[g].count().round() as usize;
        eprintln!("  {:9} {:>9} records  ~{} unique minimizers", GROUP_NAMES[g], rec_counts[g].load(Ordering::Relaxed), raw_counts[g]);
    }
    drop(hlls);
    eprintln!("  [pass A] {:.1}s", t.elapsed().as_secs_f64());

    eprintln!("\nPass B: building raw per-group filters @ {}% FPR...", RAW_FPR * 100.0);
    let t = Instant::now();
    let mut raw: Vec<Option<AtomicBloom>> = (0..NGROUPS)
        .map(|g| {
            if raw_counts[g] == 0 {
                None
            } else {
                eprintln!(" {}:", GROUP_NAMES[g]);
                Some(AtomicBloom::new(raw_counts[g], RAW_FPR))
            }
        })
        .collect();
    {
        let raw_ref = &raw;
        let dummy: Vec<AtomicU64> = (0..NGROUPS).map(|_| AtomicU64::new(0)).collect();
        stream_reference(&ref_files, &group_of, &dummy, |batch| {
            let mut mins = Vec::new();
            for (g, seq) in batch {
                if let Some(b) = &raw_ref[*g] {
                    scanner.scan_into(seq, &mut mins);
                    for &m in &mins {
                        if m != u64::MAX {
                            b.insert(m);
                        }
                    }
                }
            }
            Vec::new()
        })?;
    }
    eprintln!("  [pass B] {:.1}s", t.elapsed().as_secs_f64());

    // bacteria last so its raw filter, the largest, is freed before its output filter is built
    requested.sort_by_key(|&g| (g == BACTERIA, g));
    for &d in &requested {
        if d == BACTERIA {
            raw[BACTERIA] = None;
        }
        build_distilled(d, &raw, &scanner, &genome_dir, &db_prefix, config.fpr)?;
    }

    eprintln!("\nDone. It took {:.2}s", total_start.elapsed().as_secs_f64());
    Ok(())
}

fn build_distilled(
    d: usize,
    raw: &[Option<AtomicBloom>],
    scanner: &MinimizerScanner,
    genome_dir: &Path,
    db_prefix: &str,
    fpr: f64,
) -> Result<()> {
    let name = GROUP_NAMES[d];
    eprintln!("\nDomain {}:", name);
    let path = genome_dir.join(format!("{}.fna", name)).to_string_lossy().into_owned();
    let others: Vec<&AtomicBloom> = (0..NGROUPS).filter(|&g| g != d).filter_map(|g| raw[g].as_ref()).collect();
    let exclusive = |m: u64| !others.iter().any(|b| b.contains(m));

    let t = Instant::now();
    let total = AtomicU64::new(0);
    let kept = AtomicU64::new(0);
    let hll = stream_file(&path, |batch| {
        let mut h = new_hll();
        let mut mins = Vec::new();
        let (mut tot, mut kp) = (0u64, 0u64);
        for seq in batch {
            scanner.scan_into(seq, &mut mins);
            for &m in &mins {
                if m == u64::MAX {
                    continue;
                }
                tot += 1;
                if exclusive(m) {
                    kp += 1;
                    h.insert(&m);
                }
            }
        }
        total.fetch_add(tot, Ordering::Relaxed);
        kept.fetch_add(kp, Ordering::Relaxed);
        h
    })?;
    let mut hll = hll;
    let unique = hll.count().round() as usize;
    let (tot, kp) = (total.load(Ordering::Relaxed), kept.load(Ordering::Relaxed));
    eprintln!("  Pass C: {} minimizers (incl. duplicates), {} domain-exclusive ({:.4}% removed), ~{} unique exclusive  [{:.1}s]",
        tot, kp, 100.0 * (tot - kp) as f64 / tot.max(1) as f64, unique, t.elapsed().as_secs_f64());
    if unique == 0 {
        eprintln!("WARNING: no domain-exclusive minimizers for {}", name);
        return Ok(());
    }

    let t = Instant::now();
    let bloom = AtomicBloom::new(unique, fpr);
    {
        let bloom_ref = &bloom;
        stream_file(&path, |batch| {
            let mut mins = Vec::new();
            for seq in batch {
                scanner.scan_into(seq, &mut mins);
                for &m in &mins {
                    if m != u64::MAX && exclusive(m) {
                        bloom_ref.insert(m);
                    }
                }
            }
            new_hll()
        })?;
    }
    eprintln!("  Pass D: inserted [{:.1}s]", t.elapsed().as_secs_f64());

    let out_path = format!("{}.{}.bloom", db_prefix, name);
    let f = File::create(&out_path)?;
    let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, f);
    bincode::serialize_into(&mut writer, &bloom.into_filter())?;
    eprintln!("  Saved {}", out_path);
    Ok(())
}

struct GroupResolver {
    taxonomy: TaxonomyTree,
    prelim: FxHashMap<String, u32>,
}

impl GroupResolver {
    fn new(taxonomy: TaxonomyTree, prelim: FxHashMap<String, u32>) -> Self {
        Self { taxonomy, prelim }
    }

    fn group_of_taxid(&self, taxid: u32) -> usize {
        for t in self.taxonomy.lineage_path(taxid) {
            match t {
                TAXID_BACTERIA => return BACTERIA,
                TAXID_ARCHAEA => return ARCHAEA,
                TAXID_VIRAL => return VIRAL,
                TAXID_FUNGI => return FUNGI,
                _ => {}
            }
        }
        OTHER
    }

    fn group_of_header(&self, id: &[u8]) -> usize {
        let acc = extract_accession(id);
        let taxid = self.prelim.get(acc).copied().or_else(|| {
            acc.strip_prefix("kraken:taxid|").and_then(|r| r.split('|').next()).and_then(|s| s.parse().ok())
        });
        match taxid {
            Some(t) => self.group_of_taxid(t),
            None => OTHER,
        }
    }
}

fn stream_reference<F>(files: &[RefFile], groups: &Arc<GroupResolver>, rec_counts: &[AtomicU64], work: F) -> Result<Vec<Hll>>
where
    F: Fn(&[(usize, Vec<u8>)]) -> Vec<Hll> + Sync,
{
    let (tx, rx) = crossbeam_channel::bounded::<Vec<(usize, Vec<u8>)>>(32);
    let jobs: Vec<(String, Option<usize>)> = files.iter().map(|f| (f.path.clone(), f.fixed)).collect();
    let groups = Arc::clone(groups);
    let counts: Arc<Vec<AtomicU64>> = Arc::new((0..NGROUPS).map(|_| AtomicU64::new(0)).collect());
    let counts_reader = Arc::clone(&counts);
    let reader = std::thread::spawn(move || -> Result<()> {
        for (path, fixed) in jobs {
            let mut r = needletail::parse_fastx_file(&path).map_err(|e| anyhow::anyhow!("Cannot open FASTA {}: {}", path, e))?;
            let mut batch = Vec::with_capacity(1000);
            let mut bytes = 0usize;
            while let Some(rec) = r.next() {
                let rec = rec.map_err(|e| anyhow::anyhow!("FASTA error in {}: {}", path, e))?;
                let g = fixed.unwrap_or_else(|| groups.group_of_header(rec.id()));
                counts_reader[g].fetch_add(1, Ordering::Relaxed);
                let seq = rec.seq().into_owned();
                bytes += seq.len();
                batch.push((g, seq));
                if bytes >= 8 * 1024 * 1024 {
                    if tx.send(std::mem::replace(&mut batch, Vec::with_capacity(1000))).is_err() {
                        return Ok(());
                    }
                    bytes = 0;
                }
            }
            if !batch.is_empty() {
                let _ = tx.send(batch);
            }
        }
        Ok(())
    });

    let merged = rx
        .into_iter()
        .par_bridge()
        .map(|batch| work(&batch))
        .reduce(Vec::new, merge_hll_vecs);

    match reader.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("Reader thread panicked"),
    }
    for g in 0..NGROUPS {
        rec_counts[g].fetch_add(counts[g].load(Ordering::Relaxed), Ordering::Relaxed);
    }
    if merged.is_empty() {
        return Ok((0..NGROUPS).map(|_| new_hll()).collect());
    }
    Ok(merged)
}

fn merge_hll_vecs(mut a: Vec<Hll>, b: Vec<Hll>) -> Vec<Hll> {
    if a.is_empty() {
        return b;
    }
    for (x, y) in a.iter_mut().zip(b.iter()) {
        x.merge(y).unwrap();
    }
    a
}

fn stream_file<F>(path: &str, work: F) -> Result<Hll>
where
    F: Fn(&[Vec<u8>]) -> Hll + Sync,
{
    let (tx, rx) = crossbeam_channel::bounded::<Vec<Vec<u8>>>(32);
    let p = path.to_string();
    let reader = std::thread::spawn(move || -> Result<()> {
        let mut r = needletail::parse_fastx_file(&p).map_err(|e| anyhow::anyhow!("Cannot open FASTA {}: {}", p, e))?;
        let mut batch = Vec::with_capacity(1000);
        let mut bytes = 0usize;
        while let Some(rec) = r.next() {
            let rec = rec.map_err(|e| anyhow::anyhow!("FASTA error: {}", e))?;
            let seq = rec.seq().into_owned();
            bytes += seq.len();
            batch.push(seq);
            if bytes >= 8 * 1024 * 1024 {
                if tx.send(std::mem::replace(&mut batch, Vec::with_capacity(1000))).is_err() {
                    return Ok(());
                }
                bytes = 0;
            }
        }
        if !batch.is_empty() {
            let _ = tx.send(batch);
        }
        Ok(())
    });
    let h = rx
        .into_iter()
        .par_bridge()
        .map(|batch| work(&batch))
        .reduce(new_hll, |mut a, b| {
            a.merge(&b).unwrap();
            a
        });
    match reader.join() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e),
        Err(_) => anyhow::bail!("Reader thread panicked"),
    }
    Ok(h)
}

fn load_db_meta(prefix: &str) -> Result<(usize, usize, u64, u64)> {
    let mut f = File::open(format!("{}.meta", prefix)).context("Cannot open .meta file to read scanner parameters")?;
    let mut buf4 = [0u8; 4];
    let mut buf8 = [0u8; 8];

    f.read_exact(&mut buf4)?; let k = u32::from_le_bytes(buf4) as usize;
    f.read_exact(&mut buf4)?; let l = u32::from_le_bytes(buf4) as usize;
    f.read_exact(&mut buf8)?; let spaced_seed_mask = u64::from_le_bytes(buf8);
    f.read_exact(&mut buf8)?; let toggle_mask = u64::from_le_bytes(buf8);

    Ok((k, l, spaced_seed_mask, toggle_mask))
}
