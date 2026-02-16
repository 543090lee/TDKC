use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::Path;

use anyhow::{Context, Result};
use boomphf::Mphf;

#[inline]
fn compute_fingerprint(mut kmer: u64) -> u16 {
    kmer ^= kmer >> 33;
    kmer = kmer.wrapping_mul(0xff51afd7ed558ccd);
    kmer ^= kmer >> 33;
    (kmer & 0xFFFF) as u16
}

/// Maps accession names ↔ compact u32 IDs.
pub struct AccessionRegistry {
    name_to_id: HashMap<String, u32>,
    id_to_name: Vec<String>,
}

impl AccessionRegistry {
    pub fn new() -> Self {
        Self {
            name_to_id: HashMap::new(),
            id_to_name: Vec::new(),
        }
    }

    pub fn get_or_create(&mut self, name: &str) -> u32 {
        if let Some(&id) = self.name_to_id.get(name) {
            return id;
        }
        let id = self.id_to_name.len() as u32;
        self.id_to_name.push(name.to_string());
        self.name_to_id.insert(name.to_string(), id);
        id
    }

    pub fn get_name(&self, id: u32) -> &str {
        self.id_to_name
            .get(id as usize)
            .map(|s| s.as_str())
            .unwrap_or("Unknown")
    }

    pub fn len(&self) -> usize {
        self.id_to_name.len()
    }

    pub fn save(&self, path: &str) -> Result<()> {
        let mut f = File::create(path)?;
        writeln!(f, "ID\tAccession")?;
        for (i, name) in self.id_to_name.iter().enumerate() {
            writeln!(f, "{}\t{}", i, name)?;
        }
        Ok(())
    }

    pub fn load(path: &str) -> Result<Self> {
        let file = File::open(path).context("Cannot open accession registry")?;
        let reader = BufReader::new(file);
        let mut id_to_name = Vec::new();
        let mut name_to_id = HashMap::new();

        let mut lines = reader.lines();
        lines.next(); // skip header

        for line in lines {
            let line = line?;
            if let Some(tab_pos) = line.find('\t') {
                let name = line[tab_pos + 1..].to_string();
                let id = id_to_name.len() as u32;
                name_to_id.insert(name.clone(), id);
                id_to_name.push(name);
            }
        }

        Ok(Self {
            name_to_id,
            id_to_name,
        })
    }
}


/// CSR storage for variable-length accession lists per minimizer
pub struct CsrAccessions {
    /// Offset into `data` for each minimizer index. Length = num_minimizers + 1
    offsets: Vec<u32>,
    /// Concatenated accession IDs
    data: Vec<u32>,
}

impl CsrAccessions {
    /// Build from a map of MPHF index → list of accession IDs.
    pub fn build(num_elements: usize, accessions_by_idx: &HashMap<u64, Vec<u32>>) -> Self {
        let mut offsets = Vec::with_capacity(num_elements + 1);
        let mut data = Vec::new();
        let mut current_offset: u32 = 0;

        for i in 0..num_elements {
            offsets.push(current_offset);
            if let Some(accs) = accessions_by_idx.get(&(i as u64)) {
                for &acc in accs {
                    data.push(acc);
                }
                current_offset += accs.len() as u32;
            }
        }
        offsets.push(current_offset);

        Self { offsets, data }
    }

    /// Get accession IDs for a given MPHF index.
    pub fn get(&self, idx: usize) -> &[u32] {
        if idx + 1 >= self.offsets.len() {
            return &[];
        }
        let start = self.offsets[idx] as usize;
        let end = self.offsets[idx + 1] as usize;
        &self.data[start..end]
    }

    pub fn save(&self, path: &str) -> Result<()> {
        let mut f = File::create(path)?;
        let num_elements = (self.offsets.len() - 1) as u64;
        f.write_all(&num_elements.to_le_bytes())?;
        let data_len = self.data.len() as u64;
        f.write_all(&data_len.to_le_bytes())?;

        // Write offsets
        for &off in &self.offsets {
            f.write_all(&off.to_le_bytes())?;
        }
        // Write data
        for &acc in &self.data {
            f.write_all(&acc.to_le_bytes())?;
        }
        Ok(())
    }

    pub fn load(path: &str) -> Result<Self> {
        let raw = std::fs::read(path).context("Cannot read CSR accession file")?;

        let num_elements = u64::from_le_bytes(raw[0..8].try_into().unwrap()) as usize;
        let data_len = u64::from_le_bytes(raw[8..16].try_into().unwrap()) as usize;

        let offsets_bytes = &raw[16..16 + (num_elements + 1) * 4];
        let offsets: Vec<u32> = offsets_bytes
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect();

        let data_start = 16 + (num_elements + 1) * 4;
        let data_bytes = &raw[data_start..data_start + data_len * 4];
        let acc_data: Vec<u32> = data_bytes
            .chunks_exact(4)
            .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
            .collect();

        Ok(Self { offsets, data: acc_data })
    }
}

pub struct KmerDatabase {
    pub k: usize,
    pub l: usize,
    pub spaced_seed_mask: u64,
    pub toggle_mask: u64,
    num_minimizers: usize,
    mphf: Mphf<u64>,
    fingerprints: Vec<u16>,
    taxid_indices: Vec<u8>,
    index_to_taxid: Vec<u32>,
    /// Optional accession tracking
    accessions: Option<CsrAccessions>,
}

/// Query hit result for a single minimizer.
pub struct Hit {
    pub taxid_idx: u8,
    pub accessions: Vec<u32>,
    pub is_hit: bool,
}

impl KmerDatabase {
    pub fn true_taxid(&self, idx: u8) -> u32 {
        self.index_to_taxid
            .get(idx as usize)
            .copied()
            .unwrap_or(0)
    }

    pub fn k(&self) -> usize {
        self.k
    }

    pub fn l(&self) -> usize {
        self.l
    }

    /// Query a read: extract minimizers and look up each one
    pub fn query(&self, scanner: &crate::minimizer::MinimizerScanner, seq: &[u8]) -> Vec<Hit> {
        let minimizers = scanner.scan(seq);
        let has_acc = self.accessions.is_some();

        minimizers
            .into_iter()
            .map(|m| {
                match self.mphf.try_hash(&m) {
                    Some(idx_u64) => {
                        let idx = idx_u64 as usize;
                        if idx < self.num_minimizers
                            && self.fingerprints[idx] == compute_fingerprint(m)
                        {
                            let accessions = if has_acc {
                                self.accessions
                                    .as_ref()
                                    .unwrap()
                                    .get(idx)
                                    .to_vec()
                            } else {
                                Vec::new()
                            };
                            Hit {
                                taxid_idx: self.taxid_indices[idx],
                                accessions,
                                is_hit: true,
                            }
                        } else {
                            Hit {
                                taxid_idx: 0,
                                accessions: Vec::new(),
                                is_hit: false,
                            }
                        }
                    }
                    None => Hit {
                        taxid_idx: 0,
                        accessions: Vec::new(),
                        is_hit: false,
                    },
                }
            })
            .collect()
    }

    /// Save to disk 
    pub fn save(&self, prefix: &str) -> Result<()> {
        // Meta
        {
            let mut f = File::create(format!("{}.meta", prefix))?;
            f.write_all(&(self.k as u32).to_le_bytes())?;
            f.write_all(&(self.l as u32).to_le_bytes())?;
            f.write_all(&self.spaced_seed_mask.to_le_bytes())?;
            f.write_all(&self.toggle_mask.to_le_bytes())?;
            f.write_all(&(self.num_minimizers as u64).to_le_bytes())?;
            let has_acc: u8 = if self.accessions.is_some() { 1 } else { 0 };
            f.write_all(&[has_acc])?;
        }

        // MPHF (here i used bincode for serde)
        {
            let encoded = bincode::serialize(&self.mphf)?;
            std::fs::write(format!("{}.mphf", prefix), &encoded)?;
        }

        // Fingerprints
        {
            let mut f = File::create(format!("{}.fp", prefix))?;
            for &fp in &self.fingerprints {
                f.write_all(&fp.to_le_bytes())?;
            }
        }

        // TaxID indices
        {
            let mut f = File::create(format!("{}.taxid", prefix))?;
            f.write_all(&self.taxid_indices)?;
        }

        // TaxID mapping
        {
            let mut f = File::create(format!("{}.taxmap", prefix))?;
            let sz = self.index_to_taxid.len() as u64;
            f.write_all(&sz.to_le_bytes())?;
            for &t in &self.index_to_taxid {
                f.write_all(&t.to_le_bytes())?;
            }
        }

        {
            let mut f = File::create(format!("{}.taxmap.txt", prefix))?;
            writeln!(f, "Index\tActual_TaxID")?;
            for (i, &t) in self.index_to_taxid.iter().enumerate() {
                writeln!(f, "{}\t{}", i, t)?;
            }
        }

        // Optional accession CSR
        if let Some(ref acc) = self.accessions {
            acc.save(&format!("{}.accession", prefix))?;
        }
        Ok(())
    }

    /// Load from disk.
    pub fn load(prefix: &str) -> Result<Self> {
        eprintln!("Loading database from {}", prefix);
        let load_start = std::time::Instant::now();

        // Meta
        let (k, l, spaced_seed_mask, toggle_mask, num_minimizers, has_acc) = {
            let mut f = File::open(format!("{}.meta", prefix))
                .context("Cannot open .meta file")?;
            let mut buf4 = [0u8; 4];
            let mut buf8 = [0u8; 8];

            f.read_exact(&mut buf4)?;
            let k = u32::from_le_bytes(buf4) as usize;
            f.read_exact(&mut buf4)?;
            let l = u32::from_le_bytes(buf4) as usize;
            f.read_exact(&mut buf8)?;
            let spaced_seed_mask = u64::from_le_bytes(buf8);
            f.read_exact(&mut buf8)?;
            let toggle_mask = u64::from_le_bytes(buf8);
            f.read_exact(&mut buf8)?;
            let num_minimizers = u64::from_le_bytes(buf8) as usize;
            let mut buf1 = [0u8; 1];
            f.read_exact(&mut buf1)?;
            let has_acc = buf1[0] == 1;

            (k, l, spaced_seed_mask, toggle_mask, num_minimizers, has_acc)
        };

        // MPHF
        let t = std::time::Instant::now();
        let mphf: Mphf<u64> = {
            let data = std::fs::read(format!("{}.mphf", prefix))
                .context("Cannot read .mphf file")?;
            bincode::deserialize(&data)?
        };

        // Fingerprints
        let t = std::time::Instant::now();
        let fingerprints = {
            let data = std::fs::read(format!("{}.fp", prefix))
                .context("Cannot read .fp file")?;
            assert_eq!(data.len(), num_minimizers * 2, ".fp file size mismatch");
            data.chunks_exact(2)
                .map(|c| u16::from_le_bytes([c[0], c[1]]))
                .collect::<Vec<u16>>()
        };

        // TaxID indices
        let t = std::time::Instant::now();
        let taxid_indices = {
            let mut f = File::open(format!("{}.taxid", prefix))
                .context("Cannot open .taxid file")?;
            let mut v = vec![0u8; num_minimizers];
            f.read_exact(&mut v)?;
            v
        };

        // TaxID mapping
        let t = std::time::Instant::now();
        let index_to_taxid = {
            let data = std::fs::read(format!("{}.taxmap", prefix))
                .context("Cannot read .taxmap file")?;
            let sz = u64::from_le_bytes(data[0..8].try_into().unwrap()) as usize;
            assert_eq!(data.len(), 8 + sz * 4, ".taxmap file size mismatch");
            data[8..].chunks_exact(4)
                .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                .collect::<Vec<u32>>()
        };

        let t = std::time::Instant::now();
        let acc_path = format!("{}.accession", prefix);
        let accessions = if has_acc && Path::new(&acc_path).exists() {
            Some(CsrAccessions::load(&acc_path)?)
        } else {
            None
        };
        Ok(Self {
            k,
            l,
            spaced_seed_mask,
            toggle_mask,
            num_minimizers,
            mphf,
            fingerprints,
            taxid_indices,
            index_to_taxid,
            accessions,
        })
    }
}

/// Intermediate k-mer record extracted from Kraken output.
pub struct ExtractedKmer {
    pub sequence: String,
    pub taxid: u32,
    pub accession_id: Option<u32>,
}

/// Builds a KmerDatabase from extracted k-mers.
pub struct KmerDatabaseBuilder {
    k: usize,
    l: usize,
    spaced_seed_mask: u64,
    toggle_mask: u64,
    track_accessions: bool,
}

impl KmerDatabaseBuilder {
    pub fn new(
        k: usize,
        l: usize,
        spaced_seed_mask: u64,
        toggle_mask: u64,
        track_accessions: bool,
    ) -> Self {
        Self {
            k,
            l,
            spaced_seed_mask,
            toggle_mask,
            track_accessions,
        }
    }

    pub fn build(&self, kmers: &[ExtractedKmer], _num_threads: usize) -> Result<KmerDatabase> {
        use crate::minimizer::MinimizerScanner;
        use std::collections::{BTreeSet, HashSet};

        eprintln!("\nBuilding k-mer database");

        let start = std::time::Instant::now();

        eprintln!("\nCreating TaxID mapping");
        let unique_taxids: BTreeSet<u32> = kmers.iter().map(|k| k.taxid).collect();
        eprintln!("  Found {} unique TaxIDs", unique_taxids.len());

        if unique_taxids.len() > 255 {
            anyhow::bail!(
                "There are too many unique taxids ({}) for u8 storage (max 255)",
                unique_taxids.len()
            );
        }

        let mut index_to_taxid = Vec::new();
        let mut taxid_to_index: HashMap<u32, u8> = HashMap::new();
        for (i, &taxid) in unique_taxids.iter().enumerate() {
            index_to_taxid.push(taxid);
            taxid_to_index.insert(taxid, i as u8);
        }

        eprintln!("\nExtracting minimizers");
        let scanner = MinimizerScanner::new(
            self.k,
            self.l,
            self.spaced_seed_mask,
            self.toggle_mask,
        );

        // minimizer → (taxid_index, Set<accession_id>)
        struct MinimizerInfo {
            taxid_idx: u8,
            accessions: HashSet<u32>,
        }

        let mut minimizer_map: HashMap<u64, MinimizerInfo> = HashMap::new();

        for kmer in kmers {
            if let Some(m) = scanner.first_minimizer(kmer.sequence.as_bytes()) {
                let taxid_idx = taxid_to_index[&kmer.taxid];
                let entry = minimizer_map.entry(m).or_insert_with(|| MinimizerInfo {
                    taxid_idx,
                    accessions: HashSet::new(),
                });
                // Keep first taxid seen (deterministic)
                if let Some(acc_id) = kmer.accession_id {
                    entry.accessions.insert(acc_id);
                }
            }
        }

        let num_minimizers = minimizer_map.len();
        eprintln!("  Extracted {} unique minimizers", num_minimizers);

        if num_minimizers == 0 {
            anyhow::bail!("No minimizers found");
        }

        eprintln!("\nBuilding MPHF");
        let keys: Vec<u64> = minimizer_map.keys().copied().collect();
        let mphf = Mphf::new(2.0, &keys);

        eprintln!("\nPopulating arrays");
        let mut fingerprints = vec![0u16; num_minimizers];
        let mut taxid_indices = vec![0u8; num_minimizers];
        let mut accessions_by_idx: HashMap<u64, Vec<u32>> = HashMap::new();

        for (&kmer, info) in &minimizer_map {
            let idx = mphf.hash(&kmer);
            if (idx as usize) < num_minimizers {
                fingerprints[idx as usize] = compute_fingerprint(kmer);
                taxid_indices[idx as usize] = info.taxid_idx;

                if self.track_accessions && !info.accessions.is_empty() {
                    let mut acc_list: Vec<u32> = info.accessions.iter().copied().collect();
                    acc_list.sort_unstable();
                    accessions_by_idx.insert(idx, acc_list);
                }
            }
        }

        let accessions = if self.track_accessions {
            Some(CsrAccessions::build(num_minimizers, &accessions_by_idx))
        } else {
            None
        };

        let elapsed = start.elapsed();
        eprintln!("\nDone!");
        eprintln!("Took {:.2}s", elapsed.as_secs_f64());

        Ok(KmerDatabase {
            k: self.k,
            l: self.l,
            spaced_seed_mask: self.spaced_seed_mask,
            toggle_mask: self.toggle_mask,
            num_minimizers,
            mphf,
            fingerprints,
            taxid_indices,
            index_to_taxid,
            accessions,
        })
    }
}