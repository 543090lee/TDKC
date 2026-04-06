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

    pub fn children(&self) -> &HashMap<u32, Vec<u32>> {
        &self.children
    }

    pub fn parent_map(&self) -> &HashMap<u32, u32> {
        &self.parent
    }
}

/// Kraken2-style BFS taxonomy with sequential internal IDs.
/// Nodes closer to the root get smaller IDs, enabling a fast LCA algorithm:
/// repeatedly advance the deeper node (larger ID) to its parent until they meet.
pub struct BfsTaxonomy {
    /// external taxid -> internal BFS ID
    pub ext_to_int: HashMap<u32, u32>,
    /// internal BFS ID -> external taxid
    pub int_to_ext: Vec<u32>,
    /// internal BFS ID -> parent's internal BFS ID
    pub parent: Vec<u32>,
    /// is_relevant[internal_id] = true if taxid is in (targets + all descendants)
    pub is_relevant: Vec<bool>,
}

impl BfsTaxonomy {
    /// Build from a TaxonomyTree + the set of relevant taxids (targets + descendants).
    pub fn build(tree: &TaxonomyTree, relevant_taxids: &HashSet<u32>) -> Self {
        let parent_map = tree.parent_map();
        let children_map = tree.children();

        // BFS from root (taxid 1)
        let mut ext_to_int: HashMap<u32, u32> = HashMap::new();
        let mut int_to_ext: Vec<u32> = Vec::new();
        let mut parent_int: Vec<u32> = Vec::new();

        let mut queue = VecDeque::new();
        // Assign root
        let root: u32 = 1;
        ext_to_int.insert(root, 0);
        int_to_ext.push(root);
        parent_int.push(0); // root's parent is itself
        queue.push_back(root);

        while let Some(ext_id) = queue.pop_front() {
            let my_int = ext_to_int[&ext_id];
            if let Some(kids) = children_map.get(&ext_id) {
                for &child in kids {
                    if child == ext_id {
                        continue; // skip self-loops (root)
                    }
                    let child_int = int_to_ext.len() as u32;
                    ext_to_int.insert(child, child_int);
                    int_to_ext.push(child);
                    parent_int.push(my_int);
                    queue.push_back(child);
                }
            }
        }

        // Build is_relevant vec
        let mut is_relevant = vec![false; int_to_ext.len()];
        for &ext_id in relevant_taxids {
            if let Some(&int_id) = ext_to_int.get(&ext_id) {
                is_relevant[int_id as usize] = true;
            }
        }

        eprintln!(
            "  BFS taxonomy: {} nodes, {} relevant",
            int_to_ext.len(),
            relevant_taxids.len()
        );

        Self {
            ext_to_int,
            int_to_ext,
            parent: parent_int,
            is_relevant,
        }
    }

    /// Fast LCA: walk both nodes up until they meet.
    /// Because BFS IDs are assigned top-down, the deeper node always has the larger ID.
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

    /// Check if an internal ID is in the relevant set.
    #[inline]
    pub fn is_relevant(&self, int_id: u32) -> bool {
        (int_id as usize) < self.is_relevant.len() && self.is_relevant[int_id as usize]
    }

    /// Convert external taxid to internal BFS ID. Returns None if not in tree.
    #[inline]
    pub fn to_internal(&self, ext: u32) -> Option<u32> {
        self.ext_to_int.get(&ext).copied()
    }

    /// Convert internal BFS ID to external taxid.
    #[inline]
    pub fn to_external(&self, int_id: u32) -> u32 {
        self.int_to_ext[int_id as usize]
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