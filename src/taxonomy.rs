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
            // nodes.dmp format: taxid\t|\tparent_taxid\t|\trank\t|...
            let mut parts = line.split('\t');
            let taxid: u32 = match parts.next().and_then(|s| s.trim().parse().ok()) {
                Some(v) => v,
                None => continue,
            };
            // skip separator
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

    pub fn descendants(&self, taxid: u32) -> HashSet<u32> {
        let mut result = HashSet::new();
        let mut queue = VecDeque::new();
        queue.push_back(taxid);

        while let Some(current) = queue.pop_front() {
            if !result.insert(current) {
                continue;
            }
            if let Some(kids) = self.children.get(&current) {
                for &child in kids {
                    queue.push_back(child);
                }
            }
        }

        result
    }

    /// Check if any descendants (excluding self) are in the target set.
    pub fn has_target_descendant(&self, taxid: u32, targets: &HashSet<u32>) -> bool {
        let mut queue = VecDeque::new();
        let mut visited = HashSet::new();

        // Start with children, not self
        if let Some(kids) = self.children.get(&taxid) {
            for &child in kids {
                queue.push_back(child);
            }
        }

        while let Some(current) = queue.pop_front() {
            if !visited.insert(current) {
                continue;
            }
            if targets.contains(&current) {
                return true;
            }
            if let Some(kids) = self.children.get(&current) {
                for &child in kids {
                    queue.push_back(child);
                }
            }
        }

        false
    }
}

/// Manages the mapping from arbitrary taxids to target taxids,
/// handling the rollup logic.
pub struct TargetTaxIDManager {
    /// Targets that have other targets as descendants → exact match only
    exact_match: HashSet<u32>,
    /// Any descendant taxid → the target taxid it rolls up to
    descendant_to_target: HashMap<u32, u32>,
}

impl TargetTaxIDManager {
    pub fn new(targets: &HashSet<u32>, tree: &TaxonomyTree) -> Self {
        eprintln!("\nAnalyzing target taxid hierarchy...");

        let mut exact_match = HashSet::new();
        let mut descendant_to_target = HashMap::new();

        for &target in targets {
            if tree.has_target_descendant(target, targets) {
                exact_match.insert(target);
                eprintln!(
                    "  TaxID {}: has target descendants (exact match only)",
                    target
                );
            } else {
                let descendants = tree.descendants(target);
                eprintln!(
                    "  TaxID {}: no target descendants (rolling up {} descendant taxids)",
                    target,
                    descendants.len()
                );
                for desc in descendants {
                    descendant_to_target.insert(desc, target);
                }
            }
        }

        // eprintln!("\nSummary:");
        // eprintln!(
        //     "  Targets with target descendants (exact match): {}",
        //     exact_match.len()
        // );
        // eprintln!(
        //     "  Targets without target descendants (rollup): {}",
        //     targets.len() - exact_match.len()
        // );
        // eprintln!(
        //     "  Total descendant mappings: {}",
        //     descendant_to_target.len()
        // );

        Self {
            exact_match,
            descendant_to_target,
        }
    }

    /// Check if a taxid should be included and get the target taxid to use.
    /// Returns `Some(target_taxid)` if included, `None` otherwise.
    pub fn get_target(&self, taxid: u32) -> Option<u32> {
        if self.exact_match.contains(&taxid) {
            return Some(taxid);
        }
        self.descendant_to_target.get(&taxid).copied()
    }

    /// Get all taxids that are relevant (either exact match or descendant of a rollup target).
    pub fn all_relevant_taxids(&self) -> HashSet<u32> {
        let mut all: HashSet<u32> = self.exact_match.clone();
        for &desc in self.descendant_to_target.keys() {
            all.insert(desc);
        }
        all
    }
}

/// Load target taxids from a file (one per line, with header).
pub fn load_target_taxids(path: &str) -> Result<HashSet<u32>> {
    let file = File::open(path).context("Cannot open targets file")?;
    let reader = BufReader::new(file);
    let mut targets = HashSet::new();

    let mut lines = reader.lines();
    // Skip header
    lines.next();

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