use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::Command;
use std::time::Instant;
use rustc_hash::FxHashMap;
use anyhow::{Context, Result};

use crate::taxonomy::{load_target_taxids, TargetTaxIDManager, TaxonomyTree};

pub struct PrepConfig {
    pub fasta_file: String,
    pub accession2taxid_files: Vec<String>,
    pub targets_file: String,
    pub nodes_file: String,
    pub output_dir: String,
}

pub fn run_prep(config: PrepConfig) -> Result<()> {
    let total_start = Instant::now();

    std::fs::create_dir_all(&config.output_dir)
        .context("Cannot create output directory")?;
    let (mut accessions, mut acc_to_taxid) = extract_fasta_accessions(&config.fasta_file)?;
    
    let original_accession_count = accessions.len() + acc_to_taxid.len();
    eprintln!(
        "  Found {} unique accessions ({} pre-mapped from Kraken headers)", 
        original_accession_count, 
        acc_to_taxid.len()
    );

    lookup_accession2taxid(
        &config.accession2taxid_files,
        &mut accessions,
        &mut acc_to_taxid 
    )?;

    let unmapped_count = accessions.len(); 
    eprintln!("  Total Mapped: {}", acc_to_taxid.len());
    eprintln!("  Unmapped: {} ({:.1}%)", unmapped_count, unmapped_count as f64 / original_accession_count as f64 * 100.0);
    

    let unmapped_path = format!("{}/unmapped_accessions.txt", config.output_dir);
    if !accessions.is_empty() {
        let mut w = BufWriter::new(File::create(&unmapped_path)?);
        for (_parsed_acc, full_id) in &accessions {
            write!(w, "{}\n", full_id)?;
        }
        w.flush()?;
    }

    let prelim_map_path = format!("{}/prelim_map.txt", config.output_dir);
    write_prelim_map(&prelim_map_path, &acc_to_taxid)?;

    let taxonomy = TaxonomyTree::load(&config.nodes_file)?;
    let targets = load_target_taxids(&config.targets_file)?;
    let taxid_manager = TargetTaxIDManager::new(&targets, &taxonomy);
    let relevant_taxids = taxid_manager.all_relevant_taxids();

    let target_accessions: Vec<&String> = acc_to_taxid
        .iter()
        .filter(|(_, &taxid)| relevant_taxids.contains(&taxid))
        .map(|(acc, _)| acc)
        .collect();
    eprintln!(
        "  {} accessions belong to target clades",
        target_accessions.len());

    let target_acc_path = format!("{}/target_accessions.txt", config.output_dir);
    {
        let mut w = BufWriter::new(File::create(&target_acc_path)?);
        for acc in &target_accessions {
            writeln!(w, "{}\n", acc)?;
        }
        w.flush()?;
    }
    eprintln!("  Wrote {}", target_acc_path);

    let target_fasta_path = format!("{}/target.fasta", config.output_dir);

    let output_file = File::create(&target_fasta_path)
        .context("Cannot create target.fasta")?;

    let status = Command::new("seqtk")
        .arg("subseq")
        .arg(&config.fasta_file)
        .arg(&target_acc_path)
        .stdout(output_file)
        .status()
        .context("Failed to run seqtk. double check if seqtk is properly installed.")?;

    if !status.success() {
        anyhow::bail!("seqtk subseq failed with exit code: {}", status);
    }

    eprintln!(
        "\nDone. It took {:.2}s",
        total_start.elapsed().as_secs_f64()
    );
    Ok(())
}

fn extract_fasta_accessions(
    fasta_path: &str,
) -> Result<(FxHashMap<String, String>, FxHashMap<String, u32>)> {
    let mut child = Command::new("rg")
        .arg("^>")
        .arg("--no-filename")
        .arg("--no-line-number")
        .arg(fasta_path)
        .stdout(std::process::Stdio::piped())
        .spawn()
        .context("Failed to spawn ripgrep")?;

    let stdout = child.stdout.take().context("Failed to capture stdout")?;
    let mut reader = BufReader::with_capacity(4 * 1024 * 1024, stdout);

    let mut accessions_to_lookup = FxHashMap::default();
    let mut pre_mapped_taxids = FxHashMap::default();

    let mut line = String::new();
    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim();
        if !trimmed.is_empty() {
            let header = trimmed.strip_prefix('>').unwrap_or(trimmed);
            let full_id = header.split_whitespace().next().unwrap_or(header);
            let (parsed_acc, taxid_opt) = parse_accession_and_taxid(header);
            
            if !parsed_acc.is_empty() {
                if let Some(tid) = taxid_opt {
                    pre_mapped_taxids.insert(full_id.to_string(), tid);
                } else {
                    accessions_to_lookup.insert(parsed_acc.to_string(), full_id.to_string());
                }
            }
        }
        line.clear();
    }

    let status = child.wait().context("Failed to wait on ripgrep")?;
    if !status.success() && status.code() == Some(2) {
        anyhow::bail!("ripgrep failed.");
    }
    Ok((accessions_to_lookup, pre_mapped_taxids))
}

fn parse_accession_and_taxid(header: &str) -> (&str, Option<u32>) {
    let token = header.split_whitespace().next().unwrap_or(header);
    let mut taxid = None;
    if token.starts_with("kraken:taxid|") {
        let after_prefix = &token["kraken:taxid|".len()..];
        if let Some(pipe_pos) = after_prefix.find('|') {
            if let Ok(parsed_tid) = after_prefix[..pipe_pos].parse::<u32>() {
                if parsed_tid > 0 {
                    taxid = Some(parsed_tid);
                }
            }
        }
    }
    let after_pipe = token.rsplit('|').next().unwrap_or(token);
    let parsed_acc = after_pipe.split(':').next().unwrap_or(after_pipe);
    
    (parsed_acc, taxid)
}

fn lookup_accession2taxid(
    a2t_paths: &[String],
    fasta_accessions: &mut FxHashMap<String, String>,
    acc_to_taxid: &mut FxHashMap<String, u32>,
) -> Result<()> {

    let mut remaining = fasta_accessions.len();
    if remaining == 0 {
        eprintln!("  All accessions pre-mapped from headers, skipping accession2taxid stream!");
        return Ok(());
    }

    for (file_idx, path) in a2t_paths.iter().enumerate() {
        if remaining == 0 {
            break;
        }

        let file = File::open(path)
            .map_err(|e| anyhow::anyhow!("Cannot open {}: {}", path, e))?;
        let mut reader = BufReader::with_capacity(4 * 1024 * 1024, file);

        let mut line_count = 0u64;
        let mut line = String::new();
        while reader.read_line(&mut line)? > 0 {
            line_count += 1;

            if line_count == 1 && line.starts_with("accession") {
                line.clear();
                continue;
            }

            if line_count % 10_000_000 == 0 {
                eprint!(
                    "\r    {}M lines, found {}...",
                    line_count / 1_000_000,
                    acc_to_taxid.len()
                );
            }
            let bytes = line.as_bytes();
            let mut field_start = 0;
            let mut field_idx = 0;
            let mut acc_version: &str = "";
            let mut taxid: u32 = 0;

            for (i, &b) in bytes.iter().enumerate() {
                if b == b'\t' || i == bytes.len() - 1 {
                    let end = if i == bytes.len() - 1 && b != b'\t' {
                        i + 1
                    } else {
                        i
                    };
                    
                    match field_idx {
                        1 => acc_version = &line[field_start..end],
                        2 => {
                            taxid = line[field_start..end].parse().unwrap_or(0);
                            break; 
                        }
                        _ => {}
                    }
                    field_start = i + 1;
                    field_idx += 1;
                }
            }

            if taxid > 0 {
                
                if let Some((_parsed_acc, full_id)) = fasta_accessions.remove_entry(acc_version) {
                    
                    acc_to_taxid.insert(full_id, taxid);
                    remaining -= 1;
                    
                    if remaining == 0 {
                        break;
                    }
                }
            }         
            line.clear();
        }

    }

    Ok(())
}

fn write_prelim_map(path: &str, acc_to_taxid: &FxHashMap<String, u32>) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "accession\ttaxid")?;
    for (acc, taxid) in acc_to_taxid {
        writeln!(w, "{}\t{}", acc, taxid)?;
    }
    w.flush().map_err(Into::into)
}

pub fn load_prelim_map(path: &str) -> Result<HashMap<String, u32>> {
    let file = File::open(path).context("Cannot open prelim_map.txt")?;
    let reader = BufReader::with_capacity(4 * 1024 * 1024, file);
    let mut map = HashMap::new();

    for (i, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        if i == 0 && line.starts_with("accession") {
            continue;
        }
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(tab_pos) = line.find('\t') {
            let acc = &line[..tab_pos];
            let taxid: u32 = line[tab_pos + 1..].parse().unwrap_or(0);
            if taxid > 0 {
                map.insert(acc.to_string(), taxid);
            }
        }
    }
    Ok(map)
}