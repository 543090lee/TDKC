<div align="center">

<img src="./tdkc_logo.png" width="500" alt="TDKC Logo"/>

# TDKC

### Target Distilled K-mer Classifier: Ultrafast and Memory-Efficient Metagenomic Sequence Classification for Target Pathogen Diagnostics

[![Build](https://github.com/543090lee/TDKC/actions/workflows/rust.yml/badge.svg)](https://github.com/543090lee/TDKC/actions/workflows/rust.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/built%20with-Rust-orange?logo=rust)](https://www.rust-lang.org/)
[![Paper](https://img.shields.io/badge/paper-preprint-purple)](https://github.com/543090lee/TDKC)
[![version](https://img.shields.io/badge/version-v0.1.0-lightgrey?style=flat-square&labelColor=21262d&color=30363d)](https://github.com/543090lee/TDKC)


[Paper](https://github.com/543090lee/TDKC) · [Getting Started](#quick-start) · [Citation](#citation)

</div>
---

## 📢 News & Updates
* **[Under Development]** I am working on an even faster version of TDKC featuring multithreaded gzip decompression. This upcoming release is specifically optimized to handle massive clinical diagnostic reads with maximum efficiency. Stay tuned!

---
## Installation

```bash
git clone https://github.com/543090lee/TDKC.git
cd TDKC
conda env create -f environment.yml
conda activate tdkc
cargo install --path .

# or with mamba (faster)
mamba env create -f environment.yml
mamba activate tdkc
```

The binary will be at `./target/release/tdkc` (and also on your `PATH` inside the `tdkc` env).

---

## Quick Start

TDKC has three main steps: `prep` → `build` → `query`.

### 1. Prep

Download reference sequences from NCBI via **HTTPS** (no FTP or rsync), filter to your target taxa, and produce a local accession2taxid map. The output directory is the database, and every later command points at it.

```bash
tdkc prep \
  --domains bacteria,viral,human,archaea \
  -t data/targets.txt \
  --db tdkc_db/
```

> **Tip:** When making your own target list, it's recommended to put genus-level target taxIDs along with species-level. Many k-mers get pushed up due to conserved regions.

**Adding a user-supplied FASTA (GenBank / WGS)**

You can supplement the RefSeq downloads with a local FASTA file using `--custom`.

> **Coming soon:** Support for arbitrary custom FASTA files (sequences not registered in GenBank/WGS) is under development.

```bash
tdkc prep \
  --domains bacteria,viral,human,archaea \
  -t data/targets.txt \
  --db tdkc_db/ \
  --custom /data/my_genbank_sequences.fasta
```

| Flag | Default | Description |
|------|---------|-------------|
| `--domains` | `bacteria,viral,archaea,human` | Comma-separated RefSeq domains to download. Valid values: `bacteria`, `viral`, `archaea`, `human`, `fungi`, `invertebrate`, `plant`, `plastid`, `protozoa`. UniVec_Core is always included automatically. |
| `-t` / `--targets` | — | Path to targets file (one NCBI taxid per line, any rank). A copy is saved into the db dir as `targets.txt`. |
| `-d` / `--db` | `tdkc_db` | Database output directory. This becomes the input for every subsequent command. |
| `--custom` | — | Path to a local FASTA file to include. |
| `--concurrent-downloads` | `6` | Number of genome files being streamed from NCBI in parallel. I don't recommend going over 6, NCBI server might complain... |
| `--in-flight-chunks` | `2` | Number of dust-masking jobs to run in parallel |


### 2. Build

Distill the target k-mer index. 
```bash
tdkc build \
  --db tdkc_db/ \
  -j 32
```

Add `-a` to enable per-minimizer accession tracking (TDKC-A mode).

You can find an example `targets.txt` in the `data/` directory of this repo. The target taxa list covers the most common human respiratory and enteric viruses — feel free to use it as-is or as a starting point.

| Flag | Default | Description |
|------|---------|-------------|
| `-d` / `--db` | — | Database directory (output of `prep`) |
| `-j` | all cores | Threads |
| `-k` | `35` | Window size k |
| `-l` | `31` | Minimizer length l |
| `-a` | off | Enable accession tracking |

---

### 3. Query

Query a single sample:

```bash
tdkc query \
  --db tdkc_db/ \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -j 32 \
  -o results/sample
```

Or query an entire directory of FASTQ files at once:

```bash
tdkc query \
  --db tdkc_db/ \
  -i /data/fastq_dir/ \
  -j 32 \
  -o results/
```

TDKC auto-detects reads in the directory and processes them sequentially. Output files are written into the directory specified by `-o`.

Outputs: `results/sample.output` (per-read) and `results/sample.report`.

| Flag | Default | Description |
|------|---------|-------------|
| `-d` / `--db` | — | Database directory |
| `-1` | — | R1 FASTQ (required unless `-i` is used) |
| `-2` | — | R2 FASTQ (optional, for paired-end) |
| `-i` | — | Input directory of FASTQ files (alternative to `-1`/`-2`) |
| `-g` | `2` | Min distinct minimizer hit groups |
| `-a` | off | Output per-read accession hits (requires TDKC-A db) |
| `-b` | off | Enable domain-level detection (requires built bloom filters) |

---

### (Optional) Build Domain Bloom Filters

Domain bloom filters are opt-in: pass flags for the domains you want, and `build-domain` finds the corresponding `<db>/genome/<domain>.fna` automatically. If no flags are passed, every available domain in `genome/` is built.

```bash
# Build all available domains
tdkc build-domain --db tdkc_db/ -j 32

# Or pick specific domains
tdkc build-domain \
  --db tdkc_db/ \
  --bacteria \
  --viral \
  --archaea \
  -j 32
```

Activate at query time with `-b` to classify against broad domains too.

| Flag | Default | Description |
|------|---------|-------------|
| `-d` / `--db` | — | Database directory |
| `-p` | `0.0001` | Bloom filter false positive rate (e.g. `0.001` for 0.1%) |
| `-j` | all cores | Threads |
| `--bacteria` / `--archaea` / `--viral` / `--fungi` | off | Select domains you want to include only. If none are passed, all available are built. |

> **Note:** Lower FPR values produce more accurate domain classification but require more memory. The default 0.01% FPR is recommended for most use cases — raising it significantly (e.g. above 0.1%) may introduce spurious hits...

---

## Citation

```bibtex
@article{lee2026tdkc,
  title   = {TDKC: Ultrafast and Memory-Efficient Sequence Classification
             for Target Pathogen Diagnostics},
  author  = {Lee, Seungmo and Eskin, Eleazar},
  year    = {2026},
  url     = {https://github.com/543090lee/TDKC}
}
```

---

## License

MIT © Seungmo Lee