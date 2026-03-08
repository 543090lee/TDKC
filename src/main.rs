mod build;
mod database;
mod minimizer;
mod query;
mod taxonomy;
mod utils;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "kmer-db")]
#[command(about = "Target Distilled K-mer Classifier")]
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

        #[arg(short = 'm', long, default_value_t = 2)]
        min_distinct_minimizers: usize,

        #[arg(short = 'o', long)]
        output_prefix: String,
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
            window_size,
            minimizer_size,
        } => {
            build::run_build(build::BuildConfig {
                kraken_file: kraken,
                fasta_file: fasta,
                targets_file: targets,
                nodes_file: nodes,
                db_prefix: output,
                threads,
                track_accessions: accession,
                k: window_size,
                l: minimizer_size
            })?;
        }

        Commands::Query {
            db,
            read1,
            read2,
            threads,
            accession,
            min_distinct_minimizers,
            output_prefix,
        } => {
            query::run_query(query::QueryConfig {
                db_prefix: db,
                read1_file: read1,
                read2_file: read2,
                threads,
                use_accessions: accession,
                min_distinct_minimizers,
                output_prefix: output_prefix,
            })?;
        }
    }
    Ok(())
}