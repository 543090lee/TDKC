use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::{Context, Result};


pub const TAXID_VIRAL: u32 = 10239;
pub const TAXID_BACTERIA: u32 = 2;
pub const TAXID_ARCHAEA: u32 = 2157;
pub const TAXID_FUNGI: u32 = 4751;

// 1 is root, and it goes up when k-mer hits viruses and cellular
// remmeber, viruses is non living, so LCA root makes sense
pub const TAXID_CELLULAR: u32 = 131567;
pub const TAXID_ROOT: u32 = 1;

pub struct TaxonomyTree {
    parent: HashMap<u32, u32>,
    children: HashMap<u32, Vec<u32>>,
}

impl TaxonomyTree {
    pub fn load(path: &str) -> Result<Self> {
        let file = File::open(path).context("Cannot open nodes.dmp")?;
        let reader = BufReader::new(file);

        let mut parent = HashMap::new();
        let mut children: HashMap<u32, Vec<u32>> = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            let mut parts = line.split('\t');
            let taxid: u32 = match parts.next().and_then(|s| s.trim().parse().ok()) {
                Some(v) => v,
                None => continue,
            };
            parts.next();
            let parent_taxid: u32 = match parts.next().and_then(|s| s.trim().parse().ok()) {
                Some(v) => v,
                None => continue,
            };

            parent.insert(taxid, parent_taxid);
            children.entry(parent_taxid).or_default().push(taxid);
        }
        Ok(Self { parent, children })
    }

    pub fn lineage_path(&self, mut taxid: u32) -> Vec<u32> {
        let mut path = Vec::new();
        loop {
            path.push(taxid);
            match self.parent.get(&taxid) {
                Some(&p) if p != taxid => taxid = p,
                _ => break,
            }
        }
        path.reverse();
        path
    }

    pub fn children(&self) -> &HashMap<u32, Vec<u32>> {
        &self.children
    }
}

// This is Kraken2 style assigning internal IDs to the nodes. This is only used for tracking in full-taxon setting.
// And only when building, not querying!

pub struct BfsTaxonomy {
    pub ext_to_int: HashMap<u32, u32>,
    pub int_to_ext: Vec<u32>,
    pub parent: Vec<u32>,
    pub is_relevant: Vec<bool>,
}

impl BfsTaxonomy {
    pub fn build(tree: &TaxonomyTree, relevant_taxids: &HashSet<u32>) -> Self {
        let children_map = tree.children();

        let mut ext_to_int: HashMap<u32, u32> = HashMap::new();
        let mut int_to_ext: Vec<u32> = Vec::new();
        let mut parent_int: Vec<u32> = Vec::new();

        let mut queue = VecDeque::new();
        let root: u32 = 1;
        ext_to_int.insert(root, 0);
        int_to_ext.push(root);
        parent_int.push(0);
        queue.push_back(root);

        while let Some(ext_id) = queue.pop_front() {
            let my_int = ext_to_int[&ext_id];
            if let Some(kids) = children_map.get(&ext_id) {
                for &child in kids {
                    if child == ext_id {
                        continue;
                    }
                    let child_int = int_to_ext.len() as u32;
                    ext_to_int.insert(child, child_int);
                    int_to_ext.push(child);
                    parent_int.push(my_int);
                    queue.push_back(child);
                }
            }
        }

        let mut is_relevant = vec![false; int_to_ext.len()];
        for &ext_id in relevant_taxids {
            if let Some(&int_id) = ext_to_int.get(&ext_id) {
                is_relevant[int_id as usize] = true;
            }
        }

        Self {
            ext_to_int,
            int_to_ext,
            parent: parent_int,
            is_relevant,
        }
    }

    // Just keep walking both nodes up until they meet.
    #[inline]
    pub fn lca(&self, mut a: u32, mut b: u32) -> u32 {
        while a != b {
            if a > b {
                a = self.parent[a as usize];
            } else {
                b = self.parent[b as usize];
            }
        }
        a
    }

    #[inline]
    pub fn is_relevant(&self, int_id: u32) -> bool {
        (int_id as usize) < self.is_relevant.len() && self.is_relevant[int_id as usize]
    }

    #[inline]
    pub fn to_internal(&self, ext: u32) -> Option<u32> {
        self.ext_to_int.get(&ext).copied()
    }

    #[inline]
    pub fn to_external(&self, int_id: u32) -> u32 {
        self.int_to_ext[int_id as usize]
    }
}

// This is for rolling up taxID. This basically prevents losing kmer information that is below the target taxa
pub struct TargetTaxIDManager {
    pub target_map: HashMap<u32, u32>,
}

impl TargetTaxIDManager {
    pub fn new(targets: &HashSet<u32>, tree: &TaxonomyTree) -> Self {
        let mut target_map = HashMap::new();

        for &target in targets {
            target_map.insert(target, target);
            let mut queue = VecDeque::new();
            if let Some(kids) = tree.children.get(&target) {
                queue.extend(kids.iter().copied());
            }

            while let Some(current) = queue.pop_front() {
                if targets.contains(&current) {
                    continue;
                }

                target_map.insert(current, target);
                if let Some(kids) = tree.children.get(&current) {
                    queue.extend(kids.iter().copied());
                }
            }
        }
        Self { target_map }
    }

    pub fn get_target(&self, taxid: u32) -> Option<u32> {
        self.target_map.get(&taxid).copied()
    }

    pub fn all_relevant_taxids(&self) -> HashSet<u32> {
        self.target_map.keys().copied().collect()
    }
}

pub fn load_target_taxids(path: &str) -> Result<HashSet<u32>> {
    let file = File::open(path).context("Cannot open targets file")?;
    let reader = BufReader::new(file);
    let mut targets = HashSet::new();
    let lines = reader.lines();

    for line in lines {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Ok(taxid) = line.parse::<u32>() {
            if taxid > 1 {
                targets.insert(taxid);
            }
        }
    }
    eprintln!("Loaded {} target taxids", targets.len());
    Ok(targets)
}

#[inline]
pub fn resolve_domain_lca<'a>(taxids: impl Iterator<Item = &'a u32>) -> u32 {
    let mut hit_viral = false;
    let mut hit_cellular = false;

    for &t in taxids {
        if t == TAXID_VIRAL || t == TAXID_ROOT { 
            hit_viral = true; 
        }
        if t == TAXID_BACTERIA || t == TAXID_ARCHAEA || t == TAXID_FUNGI || t == TAXID_CELLULAR { 
            hit_cellular = true; 
        }
    }

    if hit_viral && hit_cellular { TAXID_ROOT } 
    else if hit_cellular { TAXID_CELLULAR } 
    else if hit_viral { TAXID_VIRAL } 
    else { TAXID_ROOT }
}