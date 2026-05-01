import os
import glob
import csv
import argparse
import re
from collections import defaultdict

VIRUS_MAP = {
    "InfluenzaB": "11520",
    "SarsCoV2": "2697049",
    "Rhinovirus": "147711",
    "Enterovirus": "12059",
    "H3N2": "11320",
    "H1N1": "11320",
    "RSV": "11250"
}

REVERSE_VIRUS_MAP = {v: k for k, v in VIRUS_MAP.items()}

ALIAS_MAP = {
    "FluB": "InfluenzaB",
    "InfluenzaB": "InfluenzaB",
    "SarsCoV2": "SarsCoV2",
    "CoV2": "SarsCoV2",
    "Rhinovirus": "Rhinovirus",
    "Enterovirus": "Enterovirus",
    "H3N2": "H3N2",
    "H1N1": "H1N1",
    "RSV": "RSV"
}

def load_targets(filepath):
    targets = set()
    with open(filepath, 'r') as f:
        for line in f:
            tid = line.strip()
            if tid.isdigit():
                targets.add(tid)
    targets.update(VIRUS_MAP.values())
    return targets

def load_nodes(filepath):
    parent_map = {}
    if not filepath or not os.path.exists(filepath):
        return parent_map
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.split('|')
            if len(parts) >= 2:
                child = parts[0].strip()
                parent = parts[1].strip()
                parent_map[child] = parent
    return parent_map

def load_names(filepath):
    names = {}
    if not filepath or not os.path.exists(filepath):
        return names
    with open(filepath, 'r') as f:
        for line in f:
            parts = [p.strip() for p in line.split('|')]
            if len(parts) >= 4 and parts[3] == "scientific name":
                names[parts[0]] = parts[1]
    return names

def count_total_reads_from_kraken(filepath):
    
    if not filepath or not os.path.exists(filepath):
        return None
    count = 0
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip():
                count += 1
    return count

def parse_exact_reads(filepath):
    counts = defaultdict(int)
    if not os.path.exists(filepath):
        return counts
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if not parts or parts[0] == "": continue

            if len(parts) == 3:
                if parts[0] == "Target_TaxID": continue
                taxid_field = parts[0].strip()
                try:
                    read_count = int(parts[1].strip())
                    for tid in taxid_field.split(','):
                        counts[tid] += read_count
                except ValueError: pass

            elif len(parts) >= 6 and parts[0].strip().replace('.', '').isdigit():
                taxid = parts[4].strip()
                try:
                    exact_reads = int(parts[2].strip())
                    if exact_reads > 0:
                        counts[taxid] += exact_reads
                except ValueError: pass
    return counts

def resolve_taxid(taxid, expected_taxid, panel_roots, parent_map, memo):
    memo_key = (taxid, expected_taxid)
    if memo_key in memo:
        return memo[memo_key]

    if not parent_map:
        if expected_taxid and taxid == expected_taxid: return ("TP", expected_taxid)
        if taxid in panel_roots: return ("FP", taxid)
        return (None, None)

    curr = taxid
    depth = 0
    hit_expected = False
    first_fp_target = None

    while curr != '1' and curr != '0' and curr in parent_map and depth < 100:
        if expected_taxid and curr == expected_taxid:
            hit_expected = True
            break
        if curr in panel_roots and first_fp_target is None:
            first_fp_target = curr
        if parent_map[curr] == curr: break
        curr = parent_map[curr]
        depth += 1

    if hit_expected:
        result = ("TP", expected_taxid)
    elif first_fp_target is not None:
        result = ("FP", first_fp_target)
    else:
        result = (None, None)

    memo[memo_key] = result
    return result

def get_expected_virus(filename):
    if re.search(r'(?:^|_)NC(?:_|$)', filename, re.IGNORECASE):
        return "NegativeControl"
    filename_lower = filename.lower()
    for alias in sorted(ALIAS_MAP.keys(), key=len, reverse=True):
        if alias.lower() in filename_lower:
            return ALIAS_MAP[alias]
    return None

def get_concentration(filename):
    match = re.search(r'_(D\d+)', filename)
    if match:
        return match.group(1)
    return "N/A"

def clean_filename(filepath):
    name = os.path.basename(filepath)
    while True:
        prev_name = name
        name = re.sub(r'\.(report|kreport|txt|tsv|csv)$', '', name, flags=re.IGNORECASE)
        name = re.sub(r'[\._-](kraken2?|tdkc)$', '', name, flags=re.IGNORECASE)
        if name == prev_name: break
    return name

def clean_kraken_filename(filepath):
    name = os.path.basename(filepath)
    while True:
        prev_name = name
        name = re.sub(r'\.kraken$', '', name, flags=re.IGNORECASE)
        name = re.sub(r'[\._-](kraken2?)$', '', name, flags=re.IGNORECASE)
        if name == prev_name: break
    return name

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tdkc_dir", required=True)
    parser.add_argument("-k", "--kraken_dir", required=True)
    parser.add_argument("-p", "--panel_targets", required=True)
    parser.add_argument("-d", "--nodes_dmp", required=True)
    parser.add_argument("-n", "--names_dmp", required=False)
    parser.add_argument("-o", "--out_prefix", default="LOD_comparison")
    args = parser.parse_args()
    target_taxids = load_targets(args.panel_targets)
    parent_map = load_nodes(args.nodes_dmp)
    taxid_to_name = load_names(args.names_dmp)
    tree_memo = {}

    valid_exts = ('*.report', '*.kreport', '*.txt', '*.tsv')
    tdkc_files = []
    kraken_files = []
    for ext in valid_exts:
        tdkc_files.extend(glob.glob(os.path.join(args.tdkc_dir, ext)))
        kraken_files.extend(glob.glob(os.path.join(args.kraken_dir, ext)))

    kraken_raw_files = glob.glob(os.path.join(args.kraken_dir, '*.kraken'))
    kraken_raw_map = {clean_kraken_filename(f): f for f in kraken_raw_files}

    tdkc_map = {clean_filename(f): f for f in tdkc_files}
    kraken_map = {clean_filename(f): f for f in kraken_files}
    all_samples = set(tdkc_map.keys()).union(set(kraken_map.keys()))

    detailed_results = []
    all_observed_fp_names = set()

    for sample in sorted(all_samples):
        expected_virus = get_expected_virus(sample)
        if not expected_virus: continue
        target_taxid = VIRUS_MAP.get(expected_virus)
        concentration = get_concentration(sample)
        total_reads = count_total_reads_from_kraken(kraken_raw_map.get(sample))
        t_raw_counts = parse_exact_reads(tdkc_map.get(sample, "")) if sample in tdkc_map else {}
        k_raw_counts = parse_exact_reads(kraken_map.get(sample, "")) if sample in kraken_map else {}

        for tool_name, counts in [("TDKC", t_raw_counts), ("Kraken2", k_raw_counts)]:
            if not counts: continue

            tp_reads = 0
            fps = defaultdict(int)

            for tid, count in counts.items():
                if count == 0: continue
                status, resolved_id = resolve_taxid(tid, target_taxid, target_taxids, parent_map, tree_memo)
                if status == "TP":
                    tp_reads += count
                elif status == "FP":
                    raw_name = taxid_to_name.get(resolved_id, REVERSE_VIRUS_MAP.get(resolved_id, f"TaxID_{resolved_id}"))
                    clean_name = raw_name.replace(' ', '_').replace('/', '_')
                    fps[clean_name] += count
                    all_observed_fp_names.add(clean_name)

            total_fps = sum(fps.values())
            fp_breakdown_list = [f"{v}:{c}" for v, c in sorted(fps.items())]
            fp_breakdown_str = "(" + ", ".join(fp_breakdown_list) + ")" if fp_breakdown_list else ""

            row = {
                "Sample": sample,
                "Tool": tool_name,
                "Expected_Virus": expected_virus,
                "Concentration": concentration,
                "Target_TaxID": target_taxid if target_taxid else "N/A",
                "Total_Reads": total_reads if total_reads is not None else "N/A",
                "TP_Reads": tp_reads,
                "Total_Panel_FPs": total_fps,
                "FP_Breakdown": fp_breakdown_str
            }
            for v_name, fp_count in fps.items():
                row[f"FP_{v_name}"] = fp_count

            detailed_results.append(row)

    if not detailed_results:
        return

    detailed_txt = f"{args.out_prefix}_detailed.txt"
    with open(detailed_txt, 'w') as f:
        header_str = f"{'Sample':<32} | {'Tool':<7} | {'Expected':<15} | {'Conc':<5} | {'Total Reads':<11} | {'TP (Sens)':<9} | {'Total FP':<8} | {'FP Breakdown'}"
        f.write(header_str + "\n" + "-"*120 + "\n")
        for row in detailed_results:
            f.write(f"{row['Sample'][:31]:<32} | {row['Tool']:<7} | {row['Expected_Virus']:<15} | {row['Concentration']:<5} | {str(row['Total_Reads']):<11} | {row['TP_Reads']:<9} | {row['Total_Panel_FPs']:<8} | {row['FP_Breakdown']}\n")

    sorted_fp_cols = sorted(list(all_observed_fp_names))
    csv_headers = ["Sample", "Tool", "Expected_Virus", "Concentration", "Target_TaxID", "Total_Reads", "TP_Reads", "Total_Panel_FPs", "FP_Breakdown"] + [f"FP_{v}" for v in sorted_fp_cols]

    detailed_csv = f"{args.out_prefix}_detailed.csv"
    with open(detailed_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=csv_headers)
        writer.writeheader()
        for row in detailed_results:
            safe_row = {k: row.get(k, 0) for k in csv_headers}
            safe_row.update({"Sample": row["Sample"], "Tool": row["Tool"], "Expected_Virus": row["Expected_Virus"],
                             "Concentration": row["Concentration"], "Target_TaxID": row["Target_TaxID"],
                             "Total_Reads": row["Total_Reads"], "FP_Breakdown": row["FP_Breakdown"]})
            writer.writerow(safe_row)

    agg_stats = defaultdict(lambda: {
        "TDKC": {"TP": 0, "FP": 0, "count": 0},
        "Kraken2": {"TP": 0, "FP": 0, "count": 0},
        "total_reads_sum": 0,
        "total_reads_count": 0
    })
    seen_samples = defaultdict(set)

    for row in detailed_results:
        v, t = row["Expected_Virus"], row["Tool"]
        agg_stats[v][t]["TP"] += row["TP_Reads"]
        agg_stats[v][t]["FP"] += row["Total_Panel_FPs"]
        agg_stats[v][t]["count"] += 1
        if row["Total_Reads"] != "N/A" and row["Sample"] not in seen_samples[v]:
            agg_stats[v]["total_reads_sum"] += row["Total_Reads"]
            agg_stats[v]["total_reads_count"] += 1
            seen_samples[v].add(row["Sample"])

    agg_csv = f"{args.out_prefix}_aggregated.csv"
    with open(agg_csv, 'w', newline='') as f_csv:
        writer = csv.writer(f_csv)
        writer.writerow(["Virus", "Samples_N", "Total_Reads_Sum", "Mean_Total_Reads",
                         "TDKC_Total_TP", "TDKC_Total_FP", "Kraken2_Total_TP", "Kraken2_Total_FP",
                         "TDKC_Mean_TP", "Kraken2_Mean_TP"])
        for virus, tools in sorted(agg_stats.items()):
            n_t, n_k = tools["TDKC"]["count"], tools["Kraken2"]["count"]
            n = max(n_t, n_k)
            t_tp, t_fp = tools["TDKC"]["TP"], tools["TDKC"]["FP"]
            k_tp, k_fp = tools["Kraken2"]["TP"], tools["Kraken2"]["FP"]
            t_mean_tp = round(t_tp / n_t, 1) if n_t > 0 else 0
            k_mean_tp = round(k_tp / n_k, 1) if n_k > 0 else 0
            tr_sum = tools["total_reads_sum"]
            tr_count = tools["total_reads_count"]
            mean_tr = round(tr_sum / tr_count, 1) if tr_count > 0 else "N/A"
            writer.writerow([virus, n, tr_sum if tr_count > 0 else "N/A", mean_tr,
                             t_tp, t_fp, k_tp, k_fp, t_mean_tp, k_mean_tp])

if __name__ == "__main__":
    main()