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

    // Create output directory
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
        &mut acc_to_taxid // Populate the remaining into here
    )?;

    let unmapped_count = accessions.len(); 
    eprintln!("  Total Mapped: {}", acc_to_taxid.len());
    eprintln!("  Unmapped: {} ({:.1}%)", unmapped_count, unmapped_count as f64 / original_accession_count as f64 * 100.0);
    

    let unmapped_path = format!("{}/unmapped_accessions.txt", config.output_dir);
    if !accessions.is_empty() {
        let mut w = BufWriter::new(File::create(&unmapped_path)?);
        for (_parsed_acc, full_id) in &accessions {
            // Force Unix \n line endings
            write!(w, "{}\n", full_id)?;
        }
        w.flush()?;
        eprintln!("  Wrote unmapped accessions to {}", unmapped_path);
    }

    let prelim_map_path = format!("{}/prelim_map.txt", config.output_dir);
    write_prelim_map(&prelim_map_path, &acc_to_taxid)?;
    eprintln!("  Wrote {} entries", acc_to_taxid.len());

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
        "  {} accessions belong to target clades (out of {} mapped)",
        target_accessions.len(),
        acc_to_taxid.len()
    );

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

    let seqtk_start = Instant::now();
    let output_file = File::create(&target_fasta_path)
        .context("Cannot create target.fasta")?;

    let status = Command::new("seqtk")
        .arg("subseq")
        .arg(&config.fasta_file)
        .arg(&target_acc_path)
        .stdout(output_file)
        .status()
        .context("Failed to run seqtk. Is seqtk installed and in PATH?")?;

    if !status.success() {
        anyhow::bail!("seqtk subseq failed with exit code: {}", status);
    }

    // Get file size for reporting
    let fasta_size = std::fs::metadata(&target_fasta_path)
        .map(|m| m.len())
        .unwrap_or(0);
    eprintln!(
        "  Created {} ({:.2} GB) in {:.1}s",
        target_fasta_path,
        fasta_size as f64 / 1_073_741_824.0,
        seqtk_start.elapsed().as_secs_f64()
    );

    eprintln!(
        "\nPrep complete! Total time: {:.2}s",
        total_start.elapsed().as_secs_f64()
    );
    eprintln!("\nOutputs:");
    eprintln!("  {}", prelim_map_path);
    eprintln!("  {}", target_acc_path);
    eprintln!("  {}", target_fasta_path);
    eprintln!("\nNext, run build with:");
    eprintln!(
        "  kmer-db build -f {} --target-fasta {} --prelim-map {} -t {} -n {} -o <db_prefix>",
        config.fasta_file, target_fasta_path, prelim_map_path,
        config.targets_file, config.nodes_file
    );

    Ok(())
}
// Returns (Accessions to lookup, Already mapped TaxIDs)
fn extract_fasta_accessions(
    fasta_path: &str,
) -> Result<(FxHashMap<String, String>, FxHashMap<String, u32>)> {
    let start = Instant::now();

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
            
            // Use our new dual-parser!
            let (parsed_acc, taxid_opt) = parse_accession_and_taxid(header);
            
            if !parsed_acc.is_empty() {
                if let Some(tid) = taxid_opt {
                    // We already have the TaxID! Skip the 50GB file.
                    pre_mapped_taxids.insert(full_id.to_string(), tid);
                } else {
                    // Standard header, needs lookup
                    accessions_to_lookup.insert(parsed_acc.to_string(), full_id.to_string());
                }
            }
        }
        line.clear();
    }

    let status = child.wait().context("Failed to wait on ripgrep")?;
    if !status.success() && status.code() == Some(2) {
        anyhow::bail!("ripgrep failed with exit code 2");
    }

    eprintln!(
        "  Extracted {} headers in {:.1}s",
        accessions_to_lookup.len() + pre_mapped_taxids.len(),
        start.elapsed().as_secs_f64()
    );
    Ok((accessions_to_lookup, pre_mapped_taxids))
}

fn parse_accession_and_taxid(header: &str) -> (&str, Option<u32>) {
    // 1. Get the first space-delimited token
    let token = header.split_whitespace().next().unwrap_or(header);
    
    // 2. Extract TaxID if it's in Kraken format
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
    
    // 3. Get everything after the LAST pipe for the accession
    let after_pipe = token.rsplit('|').next().unwrap_or(token);
    
    // 4. Strip trailing sequence ranges (e.g., :1-61) if they exist
    let parsed_acc = after_pipe.split(':').next().unwrap_or(after_pipe);
    
    (parsed_acc, taxid)
}

fn lookup_accession2taxid(
    a2t_paths: &[String],
    fasta_accessions: &mut FxHashMap<String, String>,
    acc_to_taxid: &mut FxHashMap<String, u32>,
) -> Result<()> {

    let mut remaining = fasta_accessions.len();
    // If everything was in Kraken format, we are already done!
    if remaining == 0 {
        eprintln!("  All accessions pre-mapped from headers, skipping accession2taxid stream!");
        return Ok(());
    }

    for (file_idx, path) in a2t_paths.iter().enumerate() {
        if remaining == 0 {
            eprintln!(
                "  All accessions mapped, skipping remaining {} file(s)",
                a2t_paths.len() - file_idx
            );
            break;
        }

        eprintln!(
            "  Streaming file {}/{}: {} ({} remaining to map)...",
            file_idx + 1,
            a2t_paths.len(),
            path,
            remaining
        );

        let file = File::open(path)
            .map_err(|e| anyhow::anyhow!("Cannot open {}: {}", path, e))?;
        let mut reader = BufReader::with_capacity(4 * 1024 * 1024, file);

        let start = Instant::now();
        let mut line_count = 0u64;
        let found_before = acc_to_taxid.len();

        // 2. ONE single String allocation for the entire 50GB file!
        let mut line = String::new();

        // 3. Use read_line instead of .lines() to reuse the buffer
        while reader.read_line(&mut line)? > 0 {
            line_count += 1;

            // Skip header line
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

            // Format: accession\taccession.version\ttaxid\tgi
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
                            break; // Stop parsing the line once we have the taxid
                        }
                        _ => {}
                    }
                    field_start = i + 1;
                    field_idx += 1;
                }
            }

            if taxid > 0 {
                
                if let Some((_parsed_acc, full_id)) = fasta_accessions.remove_entry(acc_version) {
                    
                    // Insert using the FULL FASTA ID for seqtk and the downstream build phase!
                    acc_to_taxid.insert(full_id, taxid);
                    remaining -= 1;
                    
                    if remaining == 0 {
                        eprintln!(
                            "\r    Found all remaining accessions at line {}M, stopping early.",
                            line_count / 1_000_000
                        );
                        break;
                    }
                }
            }
            
            // 5. Clear the string buffer for the next line, but keep the underlying memory capacity!
            line.clear();
        }

        let found_in_file = acc_to_taxid.len() - found_before;
        eprintln!(
            "\r    {} lines, found {} new mappings in {:.1}s",
            line_count,
            found_in_file,
            start.elapsed().as_secs_f64()
        );
    }

    Ok(())
}

// ─── Write prelim_map.txt ───────────────────────────────────────────────────

fn write_prelim_map(path: &str, acc_to_taxid: &FxHashMap<String, u32>) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "accession\ttaxid")?;
    for (acc, taxid) in acc_to_taxid {
        writeln!(w, "{}\t{}", acc, taxid)?;
    }
    w.flush().map_err(Into::into)
}

/// Load prelim_map.txt back into a HashMap.
pub fn load_prelim_map(path: &str) -> Result<HashMap<String, u32>> {
    let file = File::open(path).context("Cannot open prelim_map.txt")?;
    let reader = BufReader::with_capacity(4 * 1024 * 1024, file);
    let mut map = HashMap::new();

    for (i, line_res) in reader.lines().enumerate() {
        let line = line_res?;
        if i == 0 && line.starts_with("accession") {
            continue; // skip header
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

    eprintln!("  Loaded {} accession->taxid mappings from prelim_map", map.len());
    Ok(map)
}