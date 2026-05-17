use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;
use anyhow::{Context, Result};
use boomphf::Mphf;
use rustc_hash::FxHashMap;
use crate::hash::compute_fingerprint;
use crate::taxonomy::resolve_domain_lca;
use crate::compression::{delta_vbyte_decode, delta_vbyte_encode, delta_vbyte_skip, vbyte_decode};

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
        self.id_to_name.get(id as usize).map(|s| s.as_str()).unwrap_or("Unknown")
    }

    pub fn len(&self) -> usize {
        self.id_to_name.len()
    }

    pub fn save(&self, path: &str) -> Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(w, "ID\tAccession")?;
        for (i, name) in self.id_to_name.iter().enumerate() {
            writeln!(w, "{}\t{}", i, name)?;
        }
        w.flush().map_err(Into::into)
    }

    pub fn load_for_query(path: &str) -> Result<Self> {
        let mut id_to_name = Vec::new();

        let mut lines = BufReader::new(
            File::open(path).context("Cannot open accession registry")?
        ).lines();
        lines.next();

        for line in lines {
            let line = line?;
            if let Some(tab_pos) = line.find('\t') {
                id_to_name.push(line[tab_pos + 1..].to_string());
            }
        }
        Ok(Self { name_to_id: HashMap::new(), id_to_name })
    }
}

const RANK_BLOCK_BITS: usize = 512;
const RANK_BLOCK_WORDS: usize = RANK_BLOCK_BITS / 64; 

pub struct EqClassAccessions {
    bitset: Vec<u64>,
    rank_cache: Vec<u32>,
    dense_ids: Vec<u32>,
    num_slots: usize,
    offsets: Vec<u32>,
    blob: Vec<u8>,
    pub num_classes: usize,
}

impl EqClassAccessions {
    pub fn new_empty(num_slots: usize) -> Self {
        let num_words = (num_slots + 63) / 64;
        let encoded_empty = delta_vbyte_encode(&[]);
        Self {
            bitset: vec![0u64; num_words],
            rank_cache: Vec::new(),
            dense_ids: Vec::new(),
            num_slots,
            offsets: vec![0u32],
            blob: encoded_empty,
            num_classes: 1,
        }
    }

    fn build_rank_cache(&mut self) {
        let num_blocks = (self.bitset.len() + RANK_BLOCK_WORDS - 1) / RANK_BLOCK_WORDS;
        self.rank_cache = Vec::with_capacity(num_blocks + 1);
        let mut cumulative: u32 = 0;
        for block in 0..num_blocks {
            self.rank_cache.push(cumulative);
            let start = block * RANK_BLOCK_WORDS;
            let end = (start + RANK_BLOCK_WORDS).min(self.bitset.len());
            for w in start..end {
                cumulative += self.bitset[w].count_ones();
            }
        }
        self.rank_cache.push(cumulative);
    }

    #[inline]
    fn rank(&self, idx: usize) -> usize {
        let word = idx / 64;
        let bit = idx % 64;
        let block = word / RANK_BLOCK_WORDS;
        let mut r = self.rank_cache[block] as usize;
        let block_start = block * RANK_BLOCK_WORDS;
        for w in block_start..word {
            r += self.bitset[w].count_ones() as usize;
        }
        if bit > 0 {
            r += (self.bitset[word] & ((1u64 << bit) - 1)).count_ones() as usize;
        }
        r
    }

    #[inline]
    fn has_bit(&self, idx: usize) -> bool {
        let word = idx / 64;
        let bit = idx % 64;
        if word >= self.bitset.len() {
            return false;
        }
        self.bitset[word] & (1u64 << bit) != 0
    }

    #[inline]
    fn set_bit(&mut self, idx: usize) {
        let word = idx / 64;
        let bit = idx % 64;
        self.bitset[word] |= 1u64 << bit;
    }

    fn matches(&self, class_id: u32, acc_list: &[u32]) -> bool {
        let mut pos = self.offsets[class_id as usize] as usize;
        let count = vbyte_decode(&self.blob, &mut pos) as usize;
        if count != acc_list.len() {
            return false;
        }
        let mut prev = 0u32;
        for &expected in acc_list {
            let delta = vbyte_decode(&self.blob, &mut pos);
            prev += delta;
            if prev != expected {
                return false;
            }
        }
        true
    }

    pub fn add_accessions(
        &mut self,
        idx: usize,
        acc_list: &[u32],
        hash_to_class: &mut FxHashMap<u64, u32>,
        // This build map is temp and will be dropped before save. This is just for fast look up
        build_map: &mut FxHashMap<u32, u32>,
    ) {
        if acc_list.is_empty() {
            return;
        }

        let mut h = crate::hash::hash_slice(acc_list);
        let class_id = loop {
            match hash_to_class.get(&h) {
                Some(&id) => {
                    if self.matches(id, acc_list) {
                        break id;
                    }
                    // Hash collision... linear probe
                    h = h.wrapping_add(0x9E3779B97F4A7C15);
                }
                None => {
                    let id = self.num_classes as u32;
                    self.offsets.push(self.blob.len() as u32);
                    self.blob.extend_from_slice(&delta_vbyte_encode(acc_list));
                    hash_to_class.insert(h, id);
                    self.num_classes += 1;
                    break id;
                }
            }
        };

        self.set_bit(idx);
        build_map.insert(idx as u32, class_id);
    }

    pub fn finalize_for_save(&mut self, build_map: FxHashMap<u32, u32>) {
        let total_set: usize = self.bitset.iter().map(|w| w.count_ones() as usize).sum();
        self.dense_ids = vec![0u32; total_set];

        let mut dense_pos = 0usize;
        for word_idx in 0..self.bitset.len() {
            let mut word = self.bitset[word_idx];
            while word != 0 {
                let bit = word.trailing_zeros() as usize;
                let idx = word_idx * 64 + bit;
                if let Some(&class_id) = build_map.get(&(idx as u32)) {
                    self.dense_ids[dense_pos] = class_id;
                }
                dense_pos += 1;
                word &= word - 1;
            }
        }

        self.build_rank_cache();
    }

    #[inline]
    pub fn get_class_id(&self, idx: usize) -> Option<u32> {
        if !self.has_bit(idx) {
            return None;
        }
        let dense_pos = self.rank(idx);
        Some(self.dense_ids[dense_pos])
    }

    // get_class_id, but also returns the dense position 
    #[inline]
    pub fn get_class_id_with_pos(&self, idx: usize) -> Option<(u32, u32)> {
        if !self.has_bit(idx) {
            return None;
        }
        let dense_pos = self.rank(idx);
        Some((self.dense_ids[dense_pos], dense_pos as u32))
    }

    #[inline]
    pub fn num_dense(&self) -> usize {
        self.dense_ids.len()
    }

    // this is class ID to actual accession IDs
    pub fn decode_class(&self, class_id: u32) -> Vec<u32> {
        let cid = class_id as usize;
        if cid == 0 || cid >= self.num_classes {
            return Vec::new();
        }
        let mut pos = self.offsets[cid] as usize;
        delta_vbyte_decode(&self.blob, &mut pos)
    }

    pub fn save(&self, path: &str) -> Result<()> {
        let mut w = BufWriter::new(File::create(path)?);
        let num_entries = self.dense_ids.len() as u64;
        let num_classes = self.num_classes as u64;
        let num_slots = self.num_slots as u64;
        w.write_all(&num_entries.to_le_bytes())?;
        w.write_all(&num_classes.to_le_bytes())?;
        w.write_all(&num_slots.to_le_bytes())?;

        let width: u8 = if num_classes <= 256 { 1 } else if num_classes <= 65536 { 2 } else { 4 };
        w.write_all(&[width])?;

        let num_words = self.bitset.len() as u64;
        w.write_all(&num_words.to_le_bytes())?;
        let bitset_bytes = unsafe {
            std::slice::from_raw_parts(
                self.bitset.as_ptr() as *const u8,
                self.bitset.len() * 8,
            )
        };
        w.write_all(bitset_bytes)?;

        match width {
            1 => w.write_all(&self.dense_ids.iter().map(|&c| c as u8).collect::<Vec<_>>())?,
            2 => {
                for &c in &self.dense_ids {
                    w.write_all(&(c as u16).to_le_bytes())?;
                }
            }
            _ => {
                let id_bytes = unsafe {
                    std::slice::from_raw_parts(
                        self.dense_ids.as_ptr() as *const u8,
                        self.dense_ids.len() * 4,
                    )
                };
                w.write_all(id_bytes)?;
            }
        }

        w.write_all(&(self.blob.len() as u64).to_le_bytes())?;
        w.write_all(&self.blob)?;
        w.flush().map_err(Into::into)
    }

    pub fn load(path: &str) -> Result<Self> {
        let mut f = BufReader::new(
            File::open(path).context("Cant open eq class accession file")?
        );
        let mut buf8 = [0u8; 8];
        let mut buf1 = [0u8; 1];

        f.read_exact(&mut buf8)?;
        let num_entries = u64::from_le_bytes(buf8) as usize;
        f.read_exact(&mut buf8)?;
        let num_classes = u64::from_le_bytes(buf8) as usize;
        f.read_exact(&mut buf8)?;
        let num_slots = u64::from_le_bytes(buf8) as usize;
        f.read_exact(&mut buf1)?;
        let width = buf1[0];

        f.read_exact(&mut buf8)?;
        let num_words = u64::from_le_bytes(buf8) as usize;
        let mut bitset = vec![0u64; num_words];
        let bitset_bytes = unsafe {
            std::slice::from_raw_parts_mut(
                bitset.as_mut_ptr() as *mut u8,
                num_words * 8,
            )
        };
        f.read_exact(bitset_bytes)?;

        let dense_ids: Vec<u32> = match width {
            1 => {
                let mut bytes = vec![0u8; num_entries];
                f.read_exact(&mut bytes)?;
                bytes.into_iter().map(|b| b as u32).collect()
            }
            2 => {
                let mut bytes = vec![0u8; num_entries * 2];
                f.read_exact(&mut bytes)?;
                bytes.chunks_exact(2)
                    .map(|c| u16::from_le_bytes([c[0], c[1]]) as u32)
                    .collect()
            }
            4 => {
                let mut raw = vec![0u32; num_entries];
                let byte_slice = unsafe {
                    std::slice::from_raw_parts_mut(
                        raw.as_mut_ptr() as *mut u8,
                        num_entries * 4,
                    )
                };
                f.read_exact(byte_slice)?;
                raw
            }
            _ => anyhow::bail!("Invalid class_id width: {}", width),
        };

        f.read_exact(&mut buf8)?;
        let blob_len = u64::from_le_bytes(buf8) as usize;
        let mut blob = vec![0u8; blob_len];
        f.read_exact(&mut blob)?;

        let mut offsets = Vec::with_capacity(num_classes);
        let mut bpos: usize = 0;
        for _ in 0..num_classes {
            offsets.push(bpos as u32);
            delta_vbyte_skip(&blob, &mut bpos);
        }

        let num_blocks = (bitset.len() + RANK_BLOCK_WORDS - 1) / RANK_BLOCK_WORDS;
        let mut rank_cache = Vec::with_capacity(num_blocks + 1);
        let mut cumulative: u32 = 0;
        for block in 0..num_blocks {
            rank_cache.push(cumulative);
            let start = block * RANK_BLOCK_WORDS;
            let end = (start + RANK_BLOCK_WORDS).min(bitset.len());
            for w in start..end {
                cumulative += bitset[w].count_ones();
            }
        }
        rank_cache.push(cumulative);

        Ok(Self {
            bitset,
            rank_cache,
            dense_ids,
            num_slots,
            offsets,
            blob,
            num_classes,
        })
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
    pub index_to_taxid: Vec<u32>,
    pub ancestor_matrix: Vec<u8>,
    accessions: Option<EqClassAccessions>,
    pub names_map: Option<rustc_hash::FxHashMap<u32, String>>,
}

impl KmerDatabase {

    // This is called by build.rs 
    pub fn new(
        k: usize,
        l: usize,
        spaced_seed_mask: u64,
        toggle_mask: u64,
        num_minimizers: usize,
        mphf: Mphf<u64>,
        fingerprints: Vec<u16>,
        taxid_indices: Vec<u8>,
        index_to_taxid: Vec<u32>,
        ancestor_matrix: Vec<u8>,
        accessions: Option<EqClassAccessions>,
        names_map: Option<rustc_hash::FxHashMap<u32, String>>,
    ) -> Self {
        Self {
            k, l, spaced_seed_mask, toggle_mask,
            num_minimizers, mphf, fingerprints, taxid_indices,
            index_to_taxid, ancestor_matrix, accessions, names_map
        }
    }

    pub fn true_taxid(&self, idx: u8) -> u32 {
        self.index_to_taxid.get(idx as usize).copied().unwrap_or(0)
    }

    pub fn k(&self) -> usize { self.k }
    pub fn l(&self) -> usize { self.l }

    #[inline(always)]
    pub fn is_ancestor(&self, ancestor_idx: u8, descendant_idx: u8) -> bool {
        let n = self.index_to_taxid.len();
        let a = ancestor_idx as usize;
        let d = descendant_idx as usize;
        if a >= n || d >= n { return false; }
        self.ancestor_matrix[a * n + d] != 0
    }

    pub fn resolve_accessions(&self, class_id: u32) -> Vec<u32> {
        match self.accessions {
            Some(ref acc) => acc.decode_class(class_id),
            None => Vec::new(),
        }
    }

    #[inline]
    pub fn accessions_ref(&self) -> Option<&EqClassAccessions> {
        self.accessions.as_ref()
    }

    #[inline(always)]
    pub fn get_taxid_index(&self, taxid: u32) -> Option<u8> {
        self.index_to_taxid.iter().position(|&t| t == taxid).map(|p| p as u8)
    }

    pub fn get_taxid_name(&self, taxid: u32) -> String {
        if taxid == u32::MAX {
            return "\t".to_string();
        }

        if let Some(names_map) = &self.names_map {
            if let Some(name) = names_map.get(&taxid) {
                return name.clone();
            }
        }
        "Unknown".to_string()
    }
    
    pub fn taxid_indices(&self) -> &[u8] {
        &self.taxid_indices
    }

    pub fn query_into(
        &self,
        scanner: &crate::minimizer::MinimizerScanner,
        seq: &[u8],
        minimizer_buf: &mut Vec<u64>,
        out: &mut Vec<Hit>,
        bg_filters: &[(DomainBloomFilter, u32)],
    ) {
        scanner.scan_into(seq, minimizer_buf);
        out.reserve(minimizer_buf.len());
        
        let acc_ref = self.accessions.as_ref();
        
        let mut last_m = u64::MAX;
        let mut last_hit = Hit { 
            taxid_idx: 0, 
            accession_class_id: None,
            accession_dense_pos: None,
            is_hit: false, 
            bg_taxid: 0 
        };

        for &m in minimizer_buf.iter() {
            if m == u64::MAX {
                out.push(Hit {
                    taxid_idx: 0,
                    accession_class_id: None,
                    accession_dense_pos: None,
                    is_hit: false,
                    bg_taxid: 0,
                });
                last_m = u64::MAX;
                continue;
            }

            if m == last_m {
                out.push(last_hit);
                continue;
            }

            last_m = m;
            let mut is_hit = false;
            let mut hit_idx = 0;

            if let Some(idx_u64) = self.mphf.try_hash(&m) {
                let idx = idx_u64 as usize;
                if idx < self.num_minimizers && self.fingerprints[idx] == compute_fingerprint(m) {
                    is_hit = true;
                    hit_idx = idx;
                }
            }

            let new_hit = if is_hit {
                let (class_id_opt, dense_pos_opt) = match acc_ref.and_then(|a| a.get_class_id_with_pos(hit_idx)) {
                    Some((cid, dp)) => (Some(cid), Some(dp)),
                    None => (None, None),
                };
                Hit {
                    taxid_idx: self.taxid_indices[hit_idx],
                    accession_class_id: class_id_opt,
                    accession_dense_pos: dense_pos_opt,
                    is_hit: true,
                    bg_taxid: 0,
                }
            } else {
                let mut hit_taxids = [0u32; 4]; 
                let mut match_count = 0;

                for (filter, t_id) in bg_filters {
                    if filter.contains(m) {
                        if match_count < 4 {
                            hit_taxids[match_count] = *t_id;
                        }
                        match_count += 1;
                    }
                }

                let bg_taxid = match match_count {
                    0 => 0,
                    1 => hit_taxids[0], 
                    _ => resolve_domain_lca(hit_taxids[..match_count.min(4)].iter()),
                };

                Hit { 
                    taxid_idx: 0, 
                    accession_class_id: None,
                    accession_dense_pos: None,
                    is_hit: false, 
                    bg_taxid 
                }
            };

            last_hit = new_hit;
            out.push(new_hit);
        }
    }

    pub fn save(&self, prefix: &str) -> Result<()> {
        // .meta
        {
            let mut w = BufWriter::new(File::create(format!("{}.meta", prefix))?);
            w.write_all(&(self.k as u32).to_le_bytes())?;
            w.write_all(&(self.l as u32).to_le_bytes())?;
            w.write_all(&self.spaced_seed_mask.to_le_bytes())?;
            w.write_all(&self.toggle_mask.to_le_bytes())?;
            w.write_all(&(self.num_minimizers as u64).to_le_bytes())?;
            w.write_all(&[self.accessions.is_some() as u8])?;
            w.flush()?;
        }
        // .mphf
        {
            let encoded = bincode::serialize(&self.mphf)?;
            std::fs::write(format!("{}.mphf", prefix), &encoded)?;
        }
        // .fp
        {
            let mut w = BufWriter::new(File::create(format!("{}.fp", prefix))?);
            let fp_bytes = unsafe {
                std::slice::from_raw_parts(
                    self.fingerprints.as_ptr() as *const u8,
                    self.fingerprints.len() * 2,
                )
            };
            w.write_all(fp_bytes)?;
            w.flush()?;
        }
        // .taxid - this is in u8 btw
        {
            let mut w = BufWriter::new(File::create(format!("{}.taxid", prefix))?);
            w.write_all(&self.taxid_indices)?;
            w.flush()?;
        }
        // .taxmap - binary
        {
            let mut w = BufWriter::new(File::create(format!("{}.taxmap", prefix))?);
            w.write_all(&(self.index_to_taxid.len() as u64).to_le_bytes())?;
            let taxmap_bytes = unsafe {
                std::slice::from_raw_parts(
                    self.index_to_taxid.as_ptr() as *const u8,
                    self.index_to_taxid.len() * 4,
                )
            };
            w.write_all(taxmap_bytes)?;
            w.write_all(&self.ancestor_matrix)?;
            w.flush()?;
        }
        
        //names.dmp is saved separately and directly in build.rs

        if let Some(ref acc) = self.accessions {
            acc.save(&format!("{}.accession", prefix))?;
        }
        Ok(())
    }

    pub fn load(prefix: &str, load_accessions: bool) -> Result<Self> {
        eprintln!("Loading database from {}", prefix);

        // .meta
        let (k, l, spaced_seed_mask, toggle_mask, num_minimizers, has_acc) = {
            let mut f = File::open(format!("{}.meta", prefix))
                .context("Cannot open .meta file")?;
            let mut buf4 = [0u8; 4];
            let mut buf8 = [0u8; 8];
            let mut buf1 = [0u8; 1];

            f.read_exact(&mut buf4)?; let k = u32::from_le_bytes(buf4) as usize;
            f.read_exact(&mut buf4)?; let l = u32::from_le_bytes(buf4) as usize;
            f.read_exact(&mut buf8)?; let spaced_seed_mask = u64::from_le_bytes(buf8);
            f.read_exact(&mut buf8)?; let toggle_mask = u64::from_le_bytes(buf8);
            f.read_exact(&mut buf8)?; let num_minimizers = u64::from_le_bytes(buf8) as usize;
            f.read_exact(&mut buf1)?; let has_acc = buf1[0] == 1;
            (k, l, spaced_seed_mask, toggle_mask, num_minimizers, has_acc)
        };

        // .mphf
       let mphf: Mphf<u64> = {
            let f = File::open(format!("{}.mphf", prefix)).context("Cannot open .mphf file")?;
            let reader = BufReader::new(f);
            bincode::deserialize_from(reader)?
        };

        // .fp
        let fingerprints = {
            let mut f = File::open(format!("{}.fp", prefix)).context("Cannot open .fp file")?;
            let mut fp = vec![0u16; num_minimizers];
            let byte_slice = unsafe {
                std::slice::from_raw_parts_mut(fp.as_mut_ptr() as *mut u8, num_minimizers * 2)
            };
            f.read_exact(byte_slice)?;
            fp
        };

        // .taxid
        let taxid_indices = {
            let mut f = File::open(format!("{}.taxid", prefix)).context("Cannot open .taxid file")?;
            let mut v = vec![0u8; num_minimizers];
            f.read_exact(&mut v)?;
            v
        };

        // .taxmap
        let (index_to_taxid, ancestor_matrix) = {
            let mut f = File::open(format!("{}.taxmap", prefix)).context("Cannot open .taxmap file")?;
            let mut buf8 = [0u8; 8];
            f.read_exact(&mut buf8)?;
            let sz = u64::from_le_bytes(buf8) as usize;
            let mut v = vec![0u32; sz];
            let byte_slice = unsafe {
                std::slice::from_raw_parts_mut(v.as_mut_ptr() as *mut u8, sz * 4)
            };
            f.read_exact(byte_slice)?;
            
            let mut anc = vec![0u8; sz * sz];
            if let Err(_) = f.read_exact(&mut anc) {
                anc.fill(0);
            }
            (v, anc)
        };

        // .names
        let names_map = {
            let path = format!("{}.names", prefix);
            if let Ok(file) = File::open(&path) {
                use std::io::BufRead;
                let reader = BufReader::new(file);
                let mut map = rustc_hash::FxHashMap::default();
                for line in reader.lines().flatten() {
                    let parts: Vec<&str> = line.split('\t').collect();
                    if parts.len() == 2 {
                        if let Ok(taxid) = parts[0].parse::<u32>() {
                            map.insert(taxid, parts[1].to_string());
                        }
                    }
                }
                Some(map)
            } else {
                None
            }
        };

        let acc_path = format!("{}.accession", prefix);
        let accessions = if load_accessions && has_acc && Path::new(&acc_path).exists() {
            Some(EqClassAccessions::load(&acc_path)?)
        } else {
            None
        };
        Ok(Self {
            k, l, spaced_seed_mask, toggle_mask, num_minimizers,
            mphf, fingerprints, taxid_indices, index_to_taxid, ancestor_matrix, accessions, names_map, 
        })
    }
}

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct DomainBloomFilter {
    pub data: Vec<u8>,
    pub num_hashes: u32,
    pub num_bits: u64,
}

impl DomainBloomFilter {
    pub fn new(data: Vec<u8>, num_hashes: u32, num_bits: u64) -> Self {
        Self { data, num_hashes, num_bits }
    }

    #[inline]
    pub fn contains(&self, item: u64) -> bool {
        let (h1, h2) = crate::hash::bloom_hash_pair(item);
        for i in 0..self.num_hashes {
            let bit_pos = (h1.wrapping_add((i as u64).wrapping_mul(h2))) % self.num_bits;
            let byte_idx = (bit_pos / 8) as usize;
            let bit_mask = 1u8 << (bit_pos % 8);
            if self.data[byte_idx] & bit_mask == 0 {
                return false;
            }
        }
        true
    }

    pub fn load(path: &str) -> Result<Self> {
        let f = File::open(path)?;
        let mut reader = BufReader::new(f);
        bincode::deserialize_from(&mut reader).map_err(Into::into)
    }
}

#[derive(Clone, Copy)]
pub struct Hit {
    pub taxid_idx: u8,
    pub accession_class_id: Option<u32>,
    pub accession_dense_pos: Option<u32>,
    pub is_hit: bool,
    pub bg_taxid: u32, 
}