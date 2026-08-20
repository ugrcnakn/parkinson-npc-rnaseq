#!/usr/bin/env bash

set -euo pipefail

OUTDIR="data/raw/fastq"

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
    echo "Processing $RUN"
    echo "======================================"

    if [ -f "$OUTDIR/${RUN}.fastq" ]; then
        echo "$RUN FASTQ already exists. Skipping."
        continue
    fi

    prefetch "$RUN" --max-size u

    fasterq-dump "$RUN" \
        --threads 8 \
        --outdir "$OUTDIR" \
        --progress

    echo "$RUN completed."
done

echo
echo "=== ALL DOWNLOADS COMPLETED ==="
