<div align="center">

<img src="./tdkc_logo.png" width="500" alt="TDKC Logo"/>

# TDKC

### Target Distilled K-mer Classifier: Ultrafast and Memory-Efficient Metagenomic Sequence Classification for Target Pathogen Diagnostics

[![Build](https://img.shields.io/badge/build-passing-brightgreen)](https://github.com/543090lee/TDKC)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/built%20with-Rust-orange?logo=rust)](https://www.rust-lang.org/)
[![Paper](https://img.shields.io/badge/paper-preprint-purple)](https://github.com/543090lee/TDKC)

[Paper](https://github.com/543090lee/TDKC) · [Getting Started](#quick-start) · [Citation](#citation)

</div>

---

## Installation

**Prerequisites:**
- [Rust ≥ 1.70](https://rustup.rs/)
- [`seqtk`](https://github.com/lh3/seqtk): `conda install -c bioconda seqtk`
- [`ripgrep`](https://github.com/BurntSushi/ripgrep): `cargo install ripgrep`


```bash
git clone https://github.com/543090lee/TDKC.git
cd TDKC
cargo build --release
```

The binary will be at `./target/release/tdkc`.

---

## Quick Start

TDKC has three main steps: `prep` → `build` → `query`.

### 1. Prep

Extract target sequences and build an accession to taxid map.

```bash
tdkc prep \
  -f /data/refseq.fna.gz \
  -x /data/nucl_gb.accession2taxid \
  -t targets.txt \
  -n nodes.dmp \
  -o prep_output/
```

Outputs: `prep_output/prelim_map.txt`, `prep_output/target.fasta`

---

### 2. Build

Distill the target k-mer index.

```bash
tdkc build \
  -f /data/refseq.fna.gz \
  --target-fasta prep_output/target.fasta \
  --prelim-map   prep_output/prelim_map.txt \
  -t /data/targets.txt \
  -n nodes.dmp \
  -o my_db \
  -j 32
```

Add `-a` to enable per-minimizer accession tracking (TDKC-A mode).
You can find `targets.txt` in the `data` directory. Target taxa list is made of most common human respiratory and enteric viruses. Feel free to use it! 

| Flag | Default | Description |
|------|---------|-------------|
| `-j` | all cores | Threads |
| `-w` | `35` | Window size k |
| `-m` | `31` | Minimizer length l |
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
| `-b` | off | Enable domain Bloom filter background labels |

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

## Input Files

| File | Description |
|------|-------------|
| `refseq.fna` | Full reference FASTA (NCBI RefSeq) |
| `targets.txt` | One NCBI taxid per line (any rank — genus, species, etc.) |
| `nodes.dmp` | NCBI taxonomy `nodes.dmp` |
| `nucl_gb.accession2taxid` `nucl_wgs.accession2taxid` | NCBI accession→taxid mapping |

**RefSeq:** https://ftp.ncbi.nlm.nih.gov/genomes/refseq/  
**Viral NT:** https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNucleotide/  
**Taxonomy:** https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/


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