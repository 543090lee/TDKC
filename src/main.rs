mod build;
mod database;
mod fasta_index;
mod minimizer;
mod query;
mod taxonomy;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "kmer-db")]
#[command(about = "Distilled k-mer classifier")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Build {
        /// kraken output file
        #[arg(short = 'k', long)]
        kraken: String,

        #[arg(short = 'f', long)]
        fasta: String,

        /// target taxids txt (one per line, with header)
        #[arg(short = 't', long)]
        targets: String,

        #[arg(short = 'n', long)]
        nodes: String,

        #[arg(short = 'o', long)]
        output: String,

        #[arg(short = 'j', long, default_value_t = num_cpus::get())]
        threads: usize,

        #[arg(long)]
        accession: bool,

    },

    Query {
        #[arg(short = 'd', long)]
        db: String,

        #[arg(short = 'r', long)]
        reads: String,

        #[arg(short = 'j', long, default_value_t = 1)]
        threads: usize,

        #[arg(long)]
        accession: bool,

        #[arg(long)]
        paired: bool,

        #[arg(short = 'c', long, default_value_t = 0.3)]
        coverage: f64,
    },
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Build {
            kraken,
            fasta,
            targets,
            nodes,
            output,
            threads,
            accession,
        } => {
            build::run_build(build::BuildConfig {
                kraken_file: kraken,
                fasta_file: fasta,
                targets_file: targets,
                nodes_file: nodes,
                db_prefix: output,
                threads,
                track_accessions: accession,
            })?;
        }

        Commands::Query {
            db,
            reads,
            threads,
            accession,
            paired,
            coverage,
        } => {
            query::run_query(query::QueryConfig {
                db_prefix: db,
                reads_file: reads,
                threads,
                use_accessions: accession,
                is_paired: paired,
                coverage_threshold: coverage, 
            })?;
        }
    }
    Ok(())
}