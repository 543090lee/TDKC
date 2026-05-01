#!/usr/bin/env python3
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--kraken", required=True)
    parser.add_argument("-t", "--tdkc", required=True)
    parser.add_argument("-x", "--taxid", required=True)
    parser.add_argument("-o", "--output", required=True)

    args = parser.parse_args()
    target = args.taxid
    tdkc_data = {}

    print("Step 1: Loading unsorted TDKC data into memory (this might take a moment for large files)...")
    with open(args.tdkc, 'r') as f_tdkc:
        for line in f_tdkc:
            parts = line.strip("\n").split('\t')
            if len(parts) < 5:
                continue
            
            read_id = parts[1].split(' ')[0]
            taxid = parts[2]
            lca = parts[4]
            tdkc_data[read_id] = (taxid, lca)

    print(f"Loaded {len(tdkc_data)} reads from TDKC.")
    print("Step 2: Scanning Kraken2 data and cross-referencing...")

    with open(args.kraken, 'r') as f_krak, open(args.output, 'w') as f_out:
        f_out.write("Discordance_Type\tRead_ID\tKraken_TaxID\tTDKC_TaxID\tKraken_LCA\tTDKC_LCA\n")

        matches_found = 0
        for line in f_krak:
            parts = line.strip("\n").split('\t')
            if len(parts) < 5:
                continue

            read_id_k = parts[1].split(' ')[0]
            taxid_k = parts[2]
            lca_k = parts[4]

            if read_id_k in tdkc_data:
                taxid_t, lca_t = tdkc_data[read_id_k]

                if taxid_t == target and taxid_k != target:
                    f_out.write(f"TDKC_Matches_Kraken_Misses\t{read_id_k}\t{taxid_k}\t{taxid_t}\t{lca_k}\t{lca_t}\n")
                    matches_found += 1

                elif taxid_k == target and taxid_t != target:
                    f_out.write(f"Kraken_Matches_TDKC_Misses\t{read_id_k}\t{taxid_k}\t{taxid_t}\t{lca_k}\t{lca_t}\n")
                    matches_found += 1

if __name__ == "__main__":
    main()