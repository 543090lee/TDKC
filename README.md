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

Download reference sequences from NCBI via **HTTPS** (no FTP or rsync), filter to your target taxa, and produce a local accession2taxid map.

```bash
tdkc prep \
  --domains bacteria,viral,human,archaea \
  -t data/targets.txt \
  -o prep_output/
```

**Adding a user-supplied FASTA (GenBank / WGS)**

You can supplement the RefSeq downloads with a local FASTA file using `--custom`. 

</details>

> **Coming soon:** Support for arbitrary custom FASTA files (sequences not registered in GenBank/WGS) is under development.

```bash
tdkc prep \
  --domains bacteria,viral,human,archaea \
  -t data/targets.txt \
  -o prep_output/ \
  --custom /data/my_genbank_sequences.fasta
```


| Flag | Default | Description |
|------|---------|-------------|
| `--domains` | `bacteria,viral,archaea,human` | Comma-separated RefSeq domains to download. Valid values: `bacteria`, `viral`, `archaea`, `human`, `fungi`, `invertebrate`, `plant`, `plastid`, `protozoa`. UniVec_Core is always included automatically. |
| `-t` / `--targets` | — | Path to targets file (one NCBI taxid per line, any rank) |
| `-o` / `--output-dir` | `prep_output` | Output directory |
| `--custom` | — | Path to a local FASTA file to include. |
| `--concurrent-downloads` | `6` | Number of genome files being streamed from NCBI in parallel. I don't recommend going over 6, NCBI server might complain... |
| `--in-flight-chunks` | `2` | Number of dust-masking jobs to run in parallel |


### 2. Build

Distill the target k-mer index.

> **Coming soon:** Instead of inputting each necessary files to build, make a single dir/database (since prep phase) that will automatically detect the files.  
> This will fix issue with needing to concatenate all domain reference sequences into one all.fna.


```bash
tdkc build \
  -f /data/all.fna \
  --target-fasta prep_output/target.fasta \
  --prelim-map   prep_output/prelim_map.txt \
  -t /data/targets.txt \
  -n nodes.dmp \
  -m names.dmp \
  -o my_db \
  -j 32
```

Add `-a` to enable per-minimizer accession tracking (TDKC-A mode).
You can find `targets.txt` in the `data` directory. Target taxa list is made of most common human respiratory and enteric viruses. Feel free to use it! 

| Flag | Default | Description |
|------|---------|-------------|
| `-j` | all cores | Threads |
| `-k` | `35` | Window size k |
| `-l` | `31` | Minimizer length l |
| `-a` | off | Enable accession tracking |

---

### 3. Query

Query a single sample:

```bash
tdkc query \
  -d my_db \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -j 32 \
  -o results/sample
```

Or query an entire directory of FASTQ files at once:

```bash
tdkc query \
  -d my_db \
  -i /data/fastq_dir/ \
  -j 32 \
  -o results/
```

TDKC will auto-detect reads in the directory and process them sequentially. Output files are written into the directory specified by `-o`.

Outputs: `results/sample.output` (per-read) and `results/sample.report`.

| Flag | Default | Description |
|------|---------|-------------|
| `-1` | — | R1 FASTQ (required unless `-i` is used) |
| `-2` | — | R2 FASTQ (optional, for paired-end) |
| `-i` | — | Input directory of FASTQ files (alternative to `-1`/`-2`) |
| `-g` | `2` | Min distinct minimizer hit groups |
| `-a` | off | Output per-read accession hits (requires TDKC-A db) |
| `-b` | off | Enable domain-level detection |

---

### (Optional) Build Domain Bloom Filters

```bash
tdkc build-domain \
  -d my_db \
  --bacteria /data/bacteria.fna \
  --viral    /data/viral.fna \
  --archaea /data/archaea.fna \
  -j 32
```

Activate at query time with `-b` to classify against broad domains too.

| Flag | Default | Description |
|------|---------|-------------|
| `-p` | `0.0001` | Bloom filter false positive rate (e.g. `0.001` for 0.1%) |
| `-j` | all cores | Threads |

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