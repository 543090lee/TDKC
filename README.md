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

**Prerequisites:** [Rust ≥ 1.70](https://rustup.rs/), [`seqtk`](https://github.com/lh3/seqtk), [`ripgrep`](https://github.com/BurntSushi/ripgrep)

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

Extract target sequences and build an accession→taxid map.

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
  -t targets.txt \
  -n nodes.dmp \
  -o my_db \
  -j 32
```

Add `-a` to enable per-minimizer accession tracking (TDKC-A mode).

| Flag | Default | Description |
|------|---------|-------------|
| `-j` | all cores | Threads |
| `-w` | `35` | Window size k |
| `-m` | `31` | Minimizer length l |
| `-a` | off | Enable accession tracking |

---

### 3. Query

```bash
tdkc query \
  -d my_db \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -j 32 \
  -o results/sample
```

Outputs: `results/sample.tsv` (per-read) and `results/sample.report` (summary).

| Flag | Default | Description |
|------|---------|-------------|
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
  -j 32
```

Activate at query time with `-b` to label non-target reads by broad domain.

---

## Input Files

| File | Description |
|------|-------------|
| `refseq.fna` | Full reference FASTA (NCBI RefSeq) |
| `targets.txt` | One NCBI taxid per line (any rank — genus, species, etc.) |
| `nodes.dmp` | NCBI taxonomy `nodes.dmp` |
| `nucl_gb.accession2taxid` | NCBI accession→taxid mapping |

**RefSeq:** https://ftp.ncbi.nlm.nih.gov/genomes/refseq/  
**Viral NT:** https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNucleotide/  
**Taxonomy:** https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/


---

## Citation

```bibtex
@article{lee2025tdkc,
  title   = {TDKC: Memory-Efficient and Fast Sequence Classification
             for Target Pathogen Diagnostics},
  author  = {Lee, Seungmo and Eskin, Eleazar},
  year    = {2025},
  url     = {https://github.com/543090lee/TDKC}
}
```

---

## License

MIT © Seungmo Lee & Eleazar Eskin 