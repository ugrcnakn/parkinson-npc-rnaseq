#!/usr/bin/env bash
set -euo pipefail

GTF="data/reference/ensembl92/Homo_sapiens.GRCh38.92.gtf"
BAM_DIR="results/alignment/bam"
OUT_DIR="results/counts"

mkdir -p "$OUT_DIR"

echo "Running featureCounts..."

featureCounts \
    -T 8 \
    -s 0 \
    -t exon \
    -g gene_id \
    -a "$GTF" \
    -o "$OUT_DIR/gene_counts.txt" \
    "$BAM_DIR"/*.sorted.bam

echo "Creating clean count matrix..."

awk '
BEGIN {
    OFS="\t"
}

!/^#/ {

    row++

    if (row == 1) {

        printf "gene_id"

        for (i=7; i<=NF; i++) {

            name=$i

            sub(/^.*\//, "", name)
            sub(/\.sorted\.bam$/, "", name)

            printf OFS name
        }

        printf "\n"

    } else {

        printf $1

        for (i=7; i<=NF; i++) {
            printf OFS $i
        }

        printf "\n"
    }
}
' "$OUT_DIR/gene_counts.txt" \
> "$OUT_DIR/count_matrix.tsv"

echo "featureCounts completed."

echo -n "Genes: "
awk 'END{print NR-1}' "$OUT_DIR/count_matrix.tsv"

echo -n "Columns: "
awk -F '\t' 'NR==1{print NF}' "$OUT_DIR/count_matrix.tsv"
