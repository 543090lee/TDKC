mod download;
mod pipeline;
use std::path::PathBuf;
use std::sync::Arc;
use anyhow::{Context, Result};
use crate::taxonomy::{load_target_taxids, TargetTaxIDManager, TaxonomyTree};
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use download::*;
use pipeline::*;

pub struct PrepConfig {
    pub domains: String,
    pub targets_file: String,
    pub db_dir: String,
    pub backend: String,
    pub concurrent_downloads: usize,
    pub in_flight_chunks: usize,
    pub custom_fasta: Option<String>,
    pub no_mask: bool,
}

fn ensure_dustmasker_available() -> Result<()> {
    match std::process::Command::new("dustmasker")
        .arg("-version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
    {
        Ok(_) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Err(anyhow::anyhow!(
            "`dustmasker` not found on PATH. Install NCBI BLAST+ (e.g. `brew install blast` \
             or `conda install -c bioconda blast`), or re-run with --no-mask to skip masking."
        )),
        Err(e) => Err(e).context("failed to invoke dustmasker"),
    }
}

pub async fn run_prep(cfg: PrepConfig) -> Result<()> {
    if !cfg.no_mask {
        ensure_dustmasker_available()?;
    }

    let out_dir = PathBuf::from(&cfg.db_dir);
    let genome_dir = out_dir.join("genome");

    tokio::fs::create_dir_all(&genome_dir)
        .await
        .context("Failed to create output genome directory")?;

    tokio::fs::copy(&cfg.targets_file, out_dir.join("targets.txt"))
        .await
        .context("Failed to copy targets file into db dir")?;

    let backend = make_backend(BackendKind::Http)?;

    let tax_paths = download_taxdump(&backend, &out_dir).await?;
    let tax_dir = out_dir.join("taxonomy");


    let nodes_path = tax_paths.nodes_dmp.to_string_lossy();
    let taxonomy = TaxonomyTree::load(nodes_path.as_ref())?;
    
    let targets = load_target_taxids(&cfg.targets_file)?;
    let taxid_manager = TargetTaxIDManager::new(&targets, &taxonomy);
    
    let relevant_taxids = Arc::new(taxid_manager.all_relevant_taxids());
    
    eprintln!("Found {} target taxids in the tree.", relevant_taxids.len());

    let mut domain_list: Vec<String> = cfg.domains.split(',')
        .map(|s| s.trim().to_lowercase())
        .filter(|s| !s.is_empty() && s != "none")
        .collect();

    if !domain_list.iter().any(|d| d == "univec" || d == "univec_core") {
        domain_list.push("univec_core".to_string());
    }

    if cfg.custom_fasta.is_some() || domain_list.iter().any(|d| d == "plasmid") {
        prefetch_acc2taxid(&backend, &tax_dir).await?;
    }

    let mut total_stats = PipelineStats::default();

    for domain in &domain_list {
        let done_marker = genome_dir.join(format!(".{}.done", domain));
        if done_marker.exists() {
            eprintln!("\nDomain '{}' already exists, skipping.", domain);
            continue;
        }

        eprintln!("\nProcessing Domain: {}", domain);

        let mut sources = Vec::new();
        let mut domain_custom_map: Option<Arc<FxHashMap<String, u32>>> = None;

        //univec is always downloaded
        if domain == "univec" || domain == "univec_core" {
            let meta = resolve_univec(domain);

            sources.push(GenomeSource {
                fetch_path: meta.fetch_path,
                gzipped: false,
                taxid: meta.taxid,
                assembly_accession: meta.filename,
            });

        } else if domain == "plasmid" {
            // plasmid has no assembly_summary.txt; use accession2taxid
            let cache_dir = genome_dir.join(".plasmid_cache");
            tokio::fs::create_dir_all(&cache_dir).await?;

            let files = list_refseq_catalog_files(&backend, "plasmid").await?;

            let mut wanted = std::collections::HashSet::new();
            for fname in files {
                let remote = format!("/genomes/refseq/plasmid/{}", fname);
                let local = cache_dir.join(fname.trim_end_matches(".gz"));
                let accs =
                    fetch_and_decompress_to(&backend, &remote, &local).await?;
                for a in accs {
                    wanted.insert(a);
                }
                sources.push(GenomeSource {
                    fetch_path: local.to_string_lossy().into_owned(),
                    gzipped: false,
                    taxid: 0, 
                    assembly_accession: fname,
                });
            }

            let map = fetch_acc2taxid_filtered(&backend, &wanted, &tax_dir).await?;
            domain_custom_map = Some(Arc::new(map));

        } else {
            // Standard RefSeq 
            let entries = fetch_assembly_summary(&backend, domain).await?;
            
            for entry in entries {
                match entry.fna_gz_path() {
                    Ok(path) => {
                        sources.push(GenomeSource {
                            fetch_path: path,
                            gzipped: true, // RefSeq is always gzipped
                            taxid: entry.taxid,
                            assembly_accession: entry.assembly_accession,
                        });
                    }
                    Err(e) => eprintln!("Warning: Skipping genome due to path error: {}", e),
                }
            }
        }

        let (shared_writers, shared_handles) = SharedWriters::new(&genome_dir, domain)?;

        let pipeline_cfg = PipelineConfig {
            domain_name: domain.clone(),
            source_label: format!("refseq:{}", domain),
            output_genome_dir: genome_dir.clone(),
            backend: Arc::clone(&backend),
            concurrent_downloads: cfg.concurrent_downloads,
            max_in_flight_chunks: cfg.in_flight_chunks,
            batch_threshold_bytes: 64 * 1024 * 1024,
            relevant_taxids: Arc::clone(&relevant_taxids),
            shared: Arc::clone(&shared_writers),
            custom_map: domain_custom_map,
            no_mask: cfg.no_mask,
            // plasmid sequences with no accession2taxid hit are dropped
            drop_unmapped: domain == "plasmid",
        };

        let stats = pipeline::run_pipeline(pipeline_cfg, sources).await?;
        finalize_shared(shared_writers, shared_handles).await?;

        if domain == "plasmid" {
            let _ = tokio::fs::remove_dir_all(genome_dir.join(".plasmid_cache")).await;
        }

        total_stats.genomes_processed += stats.genomes_processed;
        total_stats.records_processed += stats.records_processed;
        total_stats.records_in_target += stats.records_in_target;
        total_stats.chunks_masked += stats.chunks_masked;

        tokio::fs::write(&done_marker, b"")
            .await
            .with_context(|| format!("write done marker {}", done_marker.display()))?;
    }

    if let Some(custom_path) = &cfg.custom_fasta {
        let done_marker = genome_dir.join(".custom.done");
        if done_marker.exists() {
            eprintln!("\nCustom sequences already added, skipping.");
        } else {
            eprintln!("\nProcessing custom FASTA: {}", custom_path);

            let mut wanted = std::collections::HashSet::new();
            let file = std::fs::File::open(custom_path)?;
            let reader = std::io::BufReader::new(file);
            for line in std::io::BufRead::lines(reader) {
                let l = line?;
                if l.starts_with('>') {
                    if let Some(acc) = l[1..].split_whitespace().next() {
                        wanted.insert(acc.to_string());
                    }
                }
            }

            let custom_map = fetch_acc2taxid_filtered(&backend, &wanted, &tax_dir).await?;

            let custom_source = GenomeSource {
                fetch_path: custom_path.clone(),
                gzipped: custom_path.ends_with(".gz"),
                taxid: 0,
                assembly_accession: "custom_user_file".to_string(),
            };

            let (shared_writers, shared_handles) = SharedWriters::new(&genome_dir, "custom")?;

            let pipeline_cfg = PipelineConfig {
                domain_name: "custom".to_string(),
                source_label: "custom_file".to_string(),
                output_genome_dir: genome_dir.clone(),
                backend: Arc::clone(&backend),
                concurrent_downloads: 1,
                max_in_flight_chunks: cfg.in_flight_chunks,
                batch_threshold_bytes: 64 * 1024 * 1024,
                relevant_taxids: Arc::clone(&relevant_taxids),
                shared: Arc::clone(&shared_writers),
                custom_map: Some(Arc::new(custom_map)),
                no_mask: cfg.no_mask,
                drop_unmapped: false,
            };

            let stats = run_pipeline(pipeline_cfg, vec![custom_source]).await?;
            finalize_shared(shared_writers, shared_handles).await?;

            total_stats.genomes_processed += stats.genomes_processed;
            total_stats.records_in_target += stats.records_in_target;

            tokio::fs::write(&done_marker, b"")
                .await
                .with_context(|| format!("write done marker {}", done_marker.display()))?;
        }
    }

    finalize_outputs(&out_dir, &genome_dir).await?;

    eprintln!("Prep is done!");
    eprintln!("Total target accessions: {}", total_stats.records_in_target);
    Ok(())
}

pub fn load_prelim_map(path: &str) -> Result<FxHashMap<String, u32>> {
    eprintln!("Loading prelim map from {}...", path);
    let file = File::open(path).with_context(|| format!("Failed to open {}", path))?;
    let reader = BufReader::new(file);
    let mut map = FxHashMap::default();

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with("accession") {
            continue;
        }
        
        let mut parts = line.split('\t');
        if let (Some(acc), Some(taxid_str)) = (parts.next(), parts.next()) {
            if let Ok(taxid) = taxid_str.trim().parse::<u32>() {
                map.insert(acc.to_string(), taxid);
            }
        }
    }
    Ok(map)
}