#!/usr/bin/env bash
set -euo pipefail

RAW="data/raw/fastq"
OUT="results/fastqc"

mkdir -p "$OUT"

while read -r RUN; do

    echo "FastQC: $RUN"

    fastqc \
        -t 4 \
        -o "$OUT" \
        "$RAW/${RUN}.fastq"

done < runs.txt

echo "FastQC completed."
