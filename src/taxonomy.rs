use std::collections::{HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::{Context, Result};

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
            // nodes.dmp has a format like taxid\t|\tparent_taxid\t|\trank\t|
            let mut parts = line.split('\t');
            let taxid: u32 = match parts.next().and_then(|s| s.trim().parse().ok()) {
                Some(v) => v,
                None => continue,
            };

            // skipping separator
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

    /// Walk from taxid up to root, return the path root-first.
    /// Used as a hierarchical sort key for accession ordering.
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

}
    
//roll up logic implemented here
pub struct TargetTaxIDManager {
    // Maps ANY relevant taxid (exact or descendant) to its reporting target taxid
    pub target_map: HashMap<u32, u32>,
}

impl TargetTaxIDManager {
    pub fn new(targets: &HashSet<u32>, tree: &TaxonomyTree) -> Self {
        let mut target_map = HashMap::new();

        for &target in targets {
            // 1. Every target maps to itself (exact match)
            target_map.insert(target, target);

            // 2. Queue up the immediate children for BFS
            let mut queue = VecDeque::new();
            if let Some(kids) = tree.children.get(&target) {
                queue.extend(kids.iter().copied());
            }

            // 3. Traverse descendants
            while let Some(current) = queue.pop_front() {
                // If this descendant is ALSO explicitly targeted by the user,
                // we stop traversing down this branch. That child target's own BFS 
                // will handle rolling up its specific sub-clade.
                if targets.contains(&current) {
                    continue;
                }

                // Map this non-target descendant to the current target
                target_map.insert(current, target);

                // Continue down the tree
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