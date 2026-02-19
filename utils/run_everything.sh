#!/bin/bash

CLASSIFIER="/Users/seungmolee/Desktop/metagenomic/distilled_model/rust_sample_test/target/release/kmer-db"
DB_PREFIX="mydb"
INPUT_DIR="/Users/seungmolee/Desktop/metagenomic/distilled_model/kaiser2-tecan"
OUTPUT_DIR="/Users/seungmolee/Desktop/metagenomic/distilled_model/rust_sample_test/mass_result"
THREADS="14"

mkdir -p "$OUTPUT_DIR"

for fastq in "$INPUT_DIR"/*.fastq.gz; do
    filename=$(basename "$fastq")
    

    if [[ $filename =~ ([A-H][0-9]{2}).*_(R[12])_ ]]; then
        well="${BASH_REMATCH[1]}"
        read="${BASH_REMATCH[2]}"
        output_name="${well}_${read}"
    else
        output_name="${filename%.fastq.gz}"
    fi
    
    echo "Currently running: $filename -> $output_name"
    
    "$CLASSIFIER" query -d "$DB_PREFIX" -1 "$fastq" -2 "$fastq" -j "$THREADS" -a -o "$output_name" 
    
done
