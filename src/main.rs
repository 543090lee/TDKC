mod build;
mod database;
mod minimizer;
mod prep;
mod query;
mod taxonomy;
mod utils;
mod build_domain;
mod compression;
mod hash;
mod read_finder;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "tdkc")]
#[command(about = "TDKC (Target Distilled K-mer Classifier) - Ultrafast and Memory-Efficient Metagenomic Classification for Target Pathogen Diagnostics")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {

    Prep {
        #[arg(short = 'f', long)]
        fasta: String,

        #[arg(short = 'x', long, num_args = 1.., required = true)]
        accession2taxid: Vec<String>,

        #[arg(short = 't', long)]
        targets: String,

        #[arg(short = 'n', long)]
        nodes: String,

        #[arg(short = 'o', long)]
        output_dir: String,
    },

    Build {
        #[arg(short = 'f', long)]
        fasta: String,

        #[arg(long)]
        target_fasta: String,

        #[arg(long)]
        prelim_map: String,

        #[arg(short = 't', long = "target-list")]
        targets: String,

        #[arg(short = 'n', long)]
        nodes: String,

        #[arg(short = 'o', long)]
        output: String,

        #[arg(short = 'j', long, default_value_t = num_cpus::get())]
        threads: usize,

        #[arg(short = 'k', long, default_value_t = 35)]
        window_size: usize,

        #[arg(short = 'l', long, default_value_t = 31)]
        minimizer_size: usize,

        #[arg(short = 'a', long = "accession2taxid")]
        accession: bool,

         #[arg(short = 'm', long)]
        names: String,
    },

    Query {
        #[arg(short = 'd', long)]
        db: String,

        #[arg(short = '1', long, required_unless_present = "input_dir")]
        read1: Option<String>,

        #[arg(short = '2', long)]
        read2: Option<String>,

        #[arg(short = 'i', long, required_unless_present = "read1")]
        input_dir: Option<String>,

        #[arg(short = 'j', long, default_value_t = 1)]
        threads: usize,

        #[arg(short = 'a', long = "accession-tracking", conflicts_with = "no_output")]
        accession: bool,

        // Same config as Kraken2
        #[arg(short = 'g', long, default_value_t = 2)]
        minimum_hit_groups: usize,

        // This can also be output dir name, not only prefix when just single fastq file is queried.
        #[arg(short = 'o', long)]
        output_prefix: String,

        #[arg(short = 'b', long)]
        background: bool,

        #[arg(long, help = "only generates .report")]
        no_output: bool,
    },

    BuildDomain {
        #[arg(short = 'd', long)]
        db: String,

        #[arg(short = 'j', long, default_value_t = num_cpus::get())]
        threads: usize,

        // default is 0.01% fp rate, higher than this might result in spurious hits...
        #[arg(short = 'p', long, default_value_t = 0.0001)]
        fpr: f64,

        #[arg(long)]
        bacteria: Option<String>,

        #[arg(long)]
        archaea: Option<String>,

        #[arg(long)]
        viral: Option<String>,

        #[arg(long)]
        fungi: Option<String>,
    },
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Prep {
            fasta,
            accession2taxid,
            targets,
            nodes,
            output_dir,
        } => {
            prep::run_prep(prep::PrepConfig {
                fasta_file: fasta,
                accession2taxid_files: accession2taxid,
                targets_file: targets,
                nodes_file: nodes,
                output_dir,
            })?;
        }

        Commands::Build {
            fasta,
            target_fasta,
            prelim_map,
            targets,
            nodes,
            output,
            threads,
            accession,
            window_size,
            minimizer_size,
            names,
        } => {
            build::run_build(build::BuildConfig {
                fasta_file: fasta,
                target_fasta_file: target_fasta,
                prelim_map_file: prelim_map,
                targets_file: targets,
                nodes_file: nodes,
                db_prefix: output,
                threads,
                track_accessions: accession,
                k: window_size,
                l: minimizer_size,
                names_dmp_path: names,
            })?;
        }

        Commands::Query {
            db,
            read1,
            read2,
            input_dir,
            threads,
            accession,
            minimum_hit_groups,
            output_prefix,
            background,
            no_output
        } => {
            let mut samples = Vec::new();
            if let Some(dir_path) = input_dir {
                std::fs::create_dir_all(&output_prefix)?;

                samples = read_finder::discover_samples(std::path::Path::new(&dir_path))?;
                let paired = samples.iter().filter(|s| s.r2.is_some()).count();
                let single = samples.len() - paired;
                eprintln!("Detected {} paired-end reads and {} single-end reads", paired, single);
                
            } else if let Some(r1) = read1 {
                let file_name = std::path::Path::new(&r1)
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy();
                
                samples.push(read_finder::Sample {
                    name: "read".to_string(),
                    r1: std::path::PathBuf::from(r1),
                    r2: read2.map(std::path::PathBuf::from),
                });
            }

            if samples.is_empty() {
                anyhow::bail!("No valid FASTQ files to query, check your files!");
            }

            query::run_query(query::QueryConfig {
                db_prefix: db,
                samples,
                threads,
                use_accessions: accession,
                minimum_hit_groups,
                output_prefix,
                background,
                no_output
            })?;
        }

        Commands::BuildDomain {
            db,
            threads,
            fpr, 
            bacteria,
            archaea,
            viral,
            fungi,
        } => {
            build_domain::run_build_domain(build_domain::BuildDomainConfig {
                db_prefix: db,
                threads,
                fpr, 
                bacteria,
                archaea,
                viral,
                fungi,
            })?;
        }
    }
    Ok(())
}
