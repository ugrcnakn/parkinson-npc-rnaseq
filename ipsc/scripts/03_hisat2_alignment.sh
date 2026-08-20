#!/usr/bin/env bash

set -euo pipefail

INDEX="data/reference/ensembl92/hisat2_index/grch38"
SPLICESITES="data/reference/ensembl92/hisat2_aux/splicesites.tsv"
FASTQ_DIR="data/raw/fastq"
BAM_DIR="results/alignment/bam"
LOG_DIR="results/alignment/logs"

mkdir -p "$BAM_DIR" "$LOG_DIR"

RUNS=(
SRR16117935
SRR16117936
SRR16117939
SRR16117940
SRR16117943
SRR16117944
SRR16117947
SRR16117948
SRR16117951
SRR16117952
SRR16117955
SRR16117956
)

for RUN in "${RUNS[@]}"
do
    echo
    echo "======================================"
    echo "Aligning $RUN"
    echo "======================================"

    FASTQ="$FASTQ_DIR/${RUN}.fastq"
    BAM="$BAM_DIR/${RUN}.sorted.bam"
    LOG="$LOG_DIR/${RUN}.hisat2.log"

    if [ -f "$BAM" ]; then
        echo "$BAM already exists. Skipping."
        continue
    fi

    if [ ! -s "$FASTQ" ]; then
        echo "ERROR: $FASTQ missing or empty."
        exit 1
    fi

    hisat2 \
        -p 8 \
        --dta \
        --known-splicesite-infile "$SPLICESITES" \
        -x "$INDEX" \
        -U "$FASTQ" \
        2> "$LOG" \
    | samtools sort \
        -@ 4 \
        -o "$BAM"

    samtools index "$BAM"

    echo "$RUN completed."
done

echo
echo "=== ALL ALIGNMENTS COMPLETED ==="
