#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import sys

# --- CONFIGURATION ---
# Generates IDs A01, A02 ... H12
EXPECTED_SAMPLES = [f"{chr(r)}{c:02d}" for r in range(ord('A'), ord('H')+1) for c in range(1, 13)]

# ─────────────────────────────────────────────
# PARSERS
# ─────────────────────────────────────────────

def parse_kraken_report(path):
    """Parse a standard Kraken2 .kreport into a DataFrame."""
    if not os.path.exists(path):
        return pd.DataFrame(columns=["taxid", "name", "reads_clade"])
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["percent", "reads_clade", "reads_direct", "rank_code", "taxid", "name"],
            usecols=[1, 3, 4, 5]
        )
        df["name"] = df["name"].str.strip()
        return df[["taxid", "name", "reads_clade"]]
    except Exception as e:
        print(f"Error parsing Kraken report {path}: {e}")
        return pd.DataFrame(columns=["taxid", "name", "reads_clade"])


def parse_tdkc_report(path):
    """
    Parse a TDKC report into a normalised DataFrame.

    Expected format (TSV, with header):
        Target_TaxID  Read_Count  Ratio
        12022         3797983     0.854
        Ambiguous     3342        0.001

    Returns a DataFrame with columns [taxid, reads_clade] where
    non-numeric TaxIDs (e.g. 'Ambiguous') are dropped.
    """
    if not os.path.exists(path):
        return pd.DataFrame(columns=["taxid", "reads_clade"])
    try:
        df = pd.read_csv(path, sep="\t")
        # Keep only numeric taxids (drop Ambiguous, Unclassified, etc.)
        df = df[pd.to_numeric(df["Target_TaxID"], errors="coerce").notna()].copy()
        df["taxid"] = df["Target_TaxID"].astype(int)
        df["reads_clade"] = df["Read_Count"].astype(int)
        return df[["taxid", "reads_clade"]]
    except Exception as e:
        print(f"Error parsing TDKC report {path}: {e}")
        return pd.DataFrame(columns=["taxid", "reads_clade"])


# ─────────────────────────────────────────────
# FILE MAPPING
# ─────────────────────────────────────────────

def map_files_to_samples(directory, extension, verbose=False):
    """
    Scans a directory and maps Sample IDs (A01..H12) to filenames.
    Returns: dict { 'A01': 'full_path_to_file', ... }
    """
    file_map = {}
    files = [f for f in os.listdir(directory) if f.endswith(extension)]

    for sample_id in EXPECTED_SAMPLES:
        matches = [f for f in files if sample_id in f]

        if len(matches) == 1:
            file_map[sample_id] = os.path.join(directory, matches[0])
        elif len(matches) > 1:
            if verbose:
                print(f"   ⚠️  Ambiguity for {sample_id} in {os.path.basename(directory)}: {matches}. Using {matches[0]}")
            file_map[sample_id] = os.path.join(directory, matches[0])

    return file_map


def detect_kraken_extension(directory):
    """Return the extension used by Kraken2 reports in this directory (.kreport or .report)."""
    files = os.listdir(directory)
    for ext in (".kreport", ".report"):
        if any(f.endswith(ext) for f in files):
            return ext
    return ".kreport"  # fallback


def detect_tdkc_extension(directory):
    """Return the extension used by TDKC reports in this directory (.report, .kreport, .tsv, .txt)."""
    files = os.listdir(directory)
    for ext in (".report", ".kreport", ".tsv", ".txt"):
        if any(f.endswith(ext) for f in files):
            return ext
    return ".report"  # fallback


# ─────────────────────────────────────────────
# CORE ANALYSIS
# ─────────────────────────────────────────────

def get_present_taxids(path, fmt, relevant_taxids):
    """
    Return the set of relevant taxids that have >0 reads in this report.
    fmt: 'kraken' | 'tdkc'
    """
    if fmt == "kraken":
        df = parse_kraken_report(path)
    else:
        df = parse_tdkc_report(path)

    df = df[df["taxid"].isin(relevant_taxids)]
    return set(df[df["reads_clade"] > 0]["taxid"])


def get_reads_for_taxids(path, fmt, relevant_taxids):
    """
    Return a dict {taxid: reads} for all relevant taxids in this report.
    fmt: 'kraken' | 'tdkc'
    """
    if fmt == "kraken":
        df = parse_kraken_report(path)
    else:
        df = parse_tdkc_report(path)

    df = df[df["taxid"].isin(relevant_taxids)]
    return dict(zip(df["taxid"], df["reads_clade"]))


def calculate_false_positives(standard_dir, tdkc_dir, fulltaxon_dir,
                               target_tsv, verbose=False):
    """
    Main analysis function.

    Parameters
    ----------
    standard_dir   : path to Standard (ground-truth) Kraken2 .kreport directory
    tdkc_dir       : path to TDKC report directory
    fulltaxon_dir  : path to Full-Taxon Kraken2 .kreport directory
    target_tsv     : path to targets TSV (columns: taxid, name)

    Returns
    -------
    df_heatmap     : DataFrame  – FP reads per virus × database (all samples)
    df_neg_stats   : DataFrame  – FP reads per virus × database (negative samples only)
    neg_sample_ids : list       – sample IDs that were negative in Standard
    """
    # 1. Load targets — headerless single-column file (just taxids)
    targets_df = pd.read_csv(target_tsv, sep="\t", header=None, names=["taxid"])
    targets_df["taxid"] = pd.to_numeric(targets_df["taxid"], errors="coerce")
    targets_df = targets_df.dropna(subset=["taxid"])
    targets_df["taxid"] = targets_df["taxid"].astype(int)
    relevant_taxids = set(targets_df["taxid"])
    # No name column — use taxid as display label
    taxid_to_name = {tid: str(tid) for tid in relevant_taxids}

    # 2. Detect extensions — Standard and Full-Taxon are always Kraken2, TDKC is always TDKC
    std_ext  = detect_kraken_extension(standard_dir);  std_fmt  = 'kraken'
    tdkc_ext = detect_tdkc_extension(tdkc_dir);        tdkc_fmt = 'tdkc'
    ft_ext   = detect_kraken_extension(fulltaxon_dir); ft_fmt   = 'kraken'

    print(f"Standard    dir : {standard_dir}  (ext={std_ext})")
    print(f"TDKC        dir : {tdkc_dir}  (ext={tdkc_ext})")
    print(f"Full-Taxon  dir : {fulltaxon_dir}  (ext={ft_ext})")

    std_map  = map_files_to_samples(standard_dir,  std_ext,  verbose)
    tdkc_map = map_files_to_samples(tdkc_dir,      tdkc_ext, verbose)
    ft_map   = map_files_to_samples(fulltaxon_dir, ft_ext,   verbose)

    print(f"   Standard   : {len(std_map)} samples mapped")
    print(f"   TDKC       : {len(tdkc_map)} samples mapped")
    print(f"   Full-Taxon : {len(ft_map)} samples mapped")

    # 3. Accumulators — keyed by sample_id
    LABELS = ["TDKC", "Full-Taxon"]
    # fp_all/fp_neg: { sample_id: { label: total_fp_reads } }
    fp_all  = {}   # all samples
    fp_neg  = {}   # negative samples only
    neg_sample_ids = []

    test_maps = {"TDKC": (tdkc_map, tdkc_fmt), "Full-Taxon": (ft_map, ft_fmt)}

    if verbose:
        print("\n" + "-" * 60)
        print("VERBOSE MODE ON")
        print("-" * 60)

    # 4. Iterate samples
    for sample_id in EXPECTED_SAMPLES:
        if sample_id not in std_map:
            if verbose:
                print(f"⏭️  Skipping {sample_id} (missing in Standard)")
            continue

        # Ground truth
        true_taxids = get_present_taxids(std_map[sample_id], std_fmt, relevant_taxids)
        is_negative = len(true_taxids) == 0

        if is_negative:
            neg_sample_ids.append(sample_id)

        fp_all[sample_id] = {lbl: 0 for lbl in LABELS}

        if verbose:
            true_names = [taxid_to_name[t] for t in true_taxids]
            status = "NEGATIVE" if is_negative else str(true_names)
            print(f"\n📂 {sample_id}  |  Standard truth: {status}")

        # Check each test DB
        for label, (tmap, tfmt) in test_maps.items():
            if sample_id not in tmap:
                if verbose:
                    print(f"   ⚠️  {label}: missing file for {sample_id}")
                continue

            reads_dict = get_reads_for_taxids(tmap[sample_id], tfmt, relevant_taxids)

            for tid, reads in reads_dict.items():
                if reads > 0 and tid not in true_taxids:
                    fp_all[sample_id][label] += reads
                    if verbose:
                        print(f"   🔴 FP  {label}: {taxid_to_name.get(tid, tid)}  reads={reads}")
                elif reads > 0 and verbose:
                    print(f"   🟢 TP  {label}: {taxid_to_name.get(tid, tid)}  reads={reads}")

        if is_negative:
            fp_neg[sample_id] = fp_all[sample_id].copy()

    # 5. Build result DataFrames — rows=samples, cols=databases
    def build_df(fp_dict):
        rows = []
        for sid, label_counts in fp_dict.items():
            row = {"sample": sid}
            row.update(label_counts)
            rows.append(row)
        return pd.DataFrame(rows).set_index("sample")

    df_heatmap   = build_df(fp_all)
    df_neg_stats = build_df(fp_neg)

    return df_heatmap, df_neg_stats, neg_sample_ids


# ─────────────────────────────────────────────
# PLOTTING
# ─────────────────────────────────────────────

def plot_heatmap(df, output_file, title=None):
    if df is None or df.empty:
        print("No data to plot.")
        return

    df_log = np.log10(df + 1)

    num_cols = len(df.columns)
    num_rows = len(df.index)
    fig_width = max(8, num_cols * 2.5)
    fig_height = max(10, num_rows * 0.4)

    plt.figure(figsize=(fig_width, fig_height))

    ax = sns.heatmap(
        df_log,
        cmap="Reds",
        linewidths=0.5,
        annot=df,
        fmt="d",
        cbar_kws={"label": "log10(False Positive Reads + 1)"},
        annot_kws={"size": 9}
    )

    plot_title = title if title else "Total False Positive Reads"
    plt.title(plot_title, fontsize=14)
    plt.xlabel("Database", fontsize=12)
    plt.ylabel("Target Viruses", fontsize=12)
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"✅ Heatmap saved: {output_file}")


def print_neg_stats(df_neg, neg_sample_ids):
    """Pretty-print the negative-sample FP summary to stdout."""
    print("\n" + "=" * 60)
    print(f"NEGATIVE SAMPLE ANALYSIS  ({len(neg_sample_ids)} negative samples)")
    print(f"Samples: {', '.join(neg_sample_ids) if neg_sample_ids else 'none'}")
    print("=" * 60)

    if df_neg.empty:
        print("  No false positives detected in negative samples.")
        return

    total_fp = df_neg.sum(axis=0)
    print("\nTotal FP reads across all negative samples:")
    for col in df_neg.columns:
        print(f"  {col:20s}: {total_fp[col]:>10,d} reads")

    print("\nPer-sample breakdown (negative samples only):")
    print(df_neg.to_string())
    print()


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare TDKC and Full-Taxon databases against a Standard (ground-truth) "
            "Kraken2 database for false-positive reads."
        )
    )
    parser.add_argument("--standard",    "-s", required=True,
                        help="Path to Standard (ground-truth) Kraken2 .kreport directory")
    parser.add_argument("--tdkc",        "-k", required=True,
                        help="Path to TDKC report directory")
    parser.add_argument("--fulltaxon",   "-f", required=True,
                        help="Path to Full-Taxon Kraken2 .kreport directory")
    parser.add_argument("--target-tsv",  "-t", required=True,
                        help="Path to targets TSV (columns: taxid, name)")
    parser.add_argument("--output",      "-o", default="fp_heatmap.png",
                        help="Output filename for the all-samples heatmap (default: fp_heatmap.png)")
    parser.add_argument("--output-neg",  "-n", default="fp_heatmap_negatives.png",
                        help="Output filename for the negative-samples heatmap (default: fp_heatmap_negatives.png)")
    parser.add_argument("--title",       help="Custom title for the all-samples heatmap")
    parser.add_argument("--verbose",     "-v", action="store_true",
                        help="Show per-sample, per-taxon FP details")

    args = parser.parse_args()

    df_all, df_neg, neg_ids = calculate_false_positives(
        standard_dir   = args.standard,
        tdkc_dir       = args.tdkc,
        fulltaxon_dir  = args.fulltaxon,
        target_tsv     = args.target_tsv,
        verbose        = args.verbose,
    )

    # All-sample heatmap
    plot_heatmap(df_all, args.output, args.title)

    # Negative-sample stats + heatmap
    print_neg_stats(df_neg, neg_ids)
    neg_title = (args.title + " — Negative Samples Only") if args.title else "FP Reads in Negative Samples"
    plot_heatmap(df_neg, args.output_neg, neg_title)


if __name__ == "__main__":
    main()

# ─────────────────────────────────────────────
# Example usage:
#
# python plot_fp.py \
#   --standard   ../../results/Standard__twist_lod \
#   --tdkc       ../../results/TDKC__twist_lod \
#   --fulltaxon  ../../results/FullTaxon__twist_lod \
#   --target-tsv /shares/swabseq/taxa-list-clean.tsv \
#   --title "TDKC vs Full-Taxon FP Comparison"
# ─────────────────────────────────────────────