#!/bin/bash

CLASSIFIER="/Users/seungmolee/Desktop/metagenomic/distilled_model/rust_sample_test/target/release/kmer-db"
DB_PREFIX="/Users/seungmolee/Desktop/metagenomic/distilled_model/rust_sample_test/mydb"
INPUT_DIR="/Users/seungmolee/Desktop/metagenomic/distilled_model/kaiser2-tecan"
OUTPUT_DIR="/Users/seungmolee/Desktop/metagenomic/distilled_model/rust_sample_test/mass_result"
THREADS="14"

mkdir -p "$OUTPUT_DIR"

count=0

for r1 in "$INPUT_DIR"/*_R1_*.fastq.gz; do
    [ -f "$r1" ] || continue
    
    r2="${r1/_R1_/_R2_}"
    
    if [ ! -f "$r2" ]; then
        echo "No r2 found, gonna skip r2"
        continue
    fi
    
    filename=$(basename "$r1")
    if [[ $filename =~ ([A-H][0-9]{2}) ]]; then
        output_name="${BASH_REMATCH[1]}"
    else
        output_name="${filename%_R1_*.fastq.gz}"
    fi
    
    echo "Currently running $output_name"
    
    "$CLASSIFIER" query -d "$DB_PREFIX" -1 "$r1" -2 "$r2" -j "$THREADS" -a -o "$OUTPUT_DIR/${output_name}"
    
    count=$((count + 1))
done

echo "Processed $count samples."