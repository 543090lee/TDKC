mod build;
mod database;
mod minimizer;
mod prep;
mod query;
mod taxonomy;
mod utils;
mod build_domain;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "tdkc")]
#[command(about = "TDKC (Target Distilled K-mer Classifier) - Fast and Memory-Efficient Metagenomic Classification for Target Pathogen Diagnostics")]
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

        #[arg(short = 't', long)]
        targets: String,

        #[arg(short = 'n', long)]
        nodes: String,

        #[arg(short = 'o', long)]
        output: String,

        #[arg(short = 'j', long, default_value_t = num_cpus::get())]
        threads: usize,

        #[arg(short = 'w', long, default_value_t = 35)]
        window_size: usize,

        #[arg(short = 'm', long, default_value_t = 31)]
        minimizer_size: usize,

        #[arg(short = 'a', long)]
        accession: bool,
    },

    Query {
        #[arg(short = 'd', long)]
        db: String,

        #[arg(short = '1', long)]
        read1: String,

        #[arg(short = '2', long)]
        read2: Option<String>,

        #[arg(short = 'j', long, default_value_t = num_cpus::get())]
        threads: usize,

        #[arg(short = 'a', long)]
        accession: bool,

        #[arg(short = 'g', long, default_value_t = 2)]
        minimum_hit_groups: usize,

        #[arg(short = 'o', long)]
        output_prefix: String,

        #[arg(short = 'b', long)]
        background: bool,
    },

    /// Build optional probabilistic Bloom filters for background domains.
    BuildDomain {
        /// Prefix of the target database (used ONLY to read k/l parameters from the .meta file)
        #[arg(short = 'd', long)]
        db: String,

        /// Number of threads
        #[arg(short = 'j', long, default_value_t = num_cpus::get())]
        threads: usize,

        /// FASTA file for Bacteria
        #[arg(long)]
        bacteria: Option<String>,

        /// FASTA file for Archaea
        #[arg(long)]
        archaea: Option<String>,

        /// FASTA file for Viral/Plasmids
        #[arg(long)]
        viral: Option<String>,

        /// FASTA file for Fungi
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
            })?;
        }

        Commands::Query {
            db,
            read1,
            read2,
            threads,
            accession,
            minimum_hit_groups,
            output_prefix,
            background
        } => {
            query::run_query(query::QueryConfig {
                db_prefix: db,
                read1_file: read1,
                read2_file: read2,
                threads,
                use_accessions: accession,
                minimum_hit_groups,
                output_prefix,
                background
            })?;
        }

        Commands::BuildDomain {
            db,
            threads,
            bacteria,
            archaea,
            viral,
            fungi,
        } => {
            build_domain::run_build_domain(build_domain::BuildDomainConfig {
                db_prefix: db,
                threads,
                bacteria,
                archaea,
                viral,
                fungi,
            })?;
        }
    }
    Ok(())
}
