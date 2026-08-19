#!/bin/bash

PROJECT="$HOME/bioinformatics/parkinson_npc_rnaseq"

METADATA="$PROJECT/data/metadata/metadata_clean.tsv"
FASTQC_ROOT="$PROJECT/results/fastqc"

OUT="$PROJECT/results/qc_summary.tsv"
TMP_DIR="$PROJECT/tmp/qc_summary"

mkdir -p "$TMP_DIR"

# Header
printf "Run\tGroup\tReads\tLength\tGC\tDeduplicated_Percent\tPer_base_quality\tPer_sequence_quality\tPer_base_content\tGC_content\tN_content\tLength_distribution\tDuplication\tOverrepresented_sequences\tAdapter_content\n" > "$OUT"


find "$FASTQC_ROOT" -type f -name "*_fastqc.zip" | sort |
while read -r ZIP
do

    SRR=$(basename "$ZIP" _fastqc.zip)

    DATA="$TMP_DIR/${SRR}_data.txt"
    SUMMARY="$TMP_DIR/${SRR}_summary.txt"

    unzip -p "$ZIP" "${SRR}_fastqc/fastqc_data.txt" > "$DATA"
    unzip -p "$ZIP" "${SRR}_fastqc/summary.txt" > "$SUMMARY"


    # ---------------------------
    # Biological group
    # ---------------------------

    GROUP=$(awk -F'\t' -v s="$SRR" '
        NR>1 && $1==s {
            print $3
            exit
        }
    ' "$METADATA")


    # ---------------------------
    # Basic statistics
    # ---------------------------

    READS=$(awk -F'\t' '
        $1=="Total Sequences" {
            print $2
            exit
        }
    ' "$DATA")

    LENGTH=$(awk -F'\t' '
        $1=="Sequence length" {
            print $2
            exit
        }
    ' "$DATA")

    GC=$(awk -F'\t' '
        $1=="%GC" {
            print $2
            exit
        }
    ' "$DATA")

    DEDUP=$(awk '
        /^#Total Deduplicated Percentage/ {
            print $NF
            exit
        }
    ' "$DATA")


    # ---------------------------
    # FastQC module status
    # ---------------------------

    PER_BASE_QUALITY=$(awk -F'\t' '
        $2=="Per base sequence quality" {
            print $1
            exit
        }
    ' "$SUMMARY")

    PER_SEQUENCE_QUALITY=$(awk -F'\t' '
        $2=="Per sequence quality scores" {
            print $1
            exit
        }
    ' "$SUMMARY")

    PER_BASE_CONTENT=$(awk -F'\t' '
        $2=="Per base sequence content" {
            print $1
            exit
        }
    ' "$SUMMARY")

    GC_CONTENT=$(awk -F'\t' '
        $2=="Per sequence GC content" {
            print $1
            exit
        }
    ' "$SUMMARY")

    N_CONTENT=$(awk -F'\t' '
        $2=="Per base N content" {
            print $1
            exit
        }
    ' "$SUMMARY")

    LENGTH_DIST=$(awk -F'\t' '
        $2=="Sequence Length Distribution" {
            print $1
            exit
        }
    ' "$SUMMARY")

    DUPLICATION=$(awk -F'\t' '
        $2=="Sequence Duplication Levels" {
            print $1
            exit
        }
    ' "$SUMMARY")

    OVERREP=$(awk -F'\t' '
        $2=="Overrepresented sequences" {
            print $1
            exit
        }
    ' "$SUMMARY")

    ADAPTER=$(awk -F'\t' '
        $2=="Adapter Content" {
            print $1
            exit
        }
    ' "$SUMMARY")


    # ---------------------------
    # Write one row
    # ---------------------------

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$SRR" \
        "$GROUP" \
        "$READS" \
        "$LENGTH" \
        "$GC" \
        "$DEDUP" \
        "$PER_BASE_QUALITY" \
        "$PER_SEQUENCE_QUALITY" \
        "$PER_BASE_CONTENT" \
        "$GC_CONTENT" \
        "$N_CONTENT" \
        "$LENGTH_DIST" \
        "$DUPLICATION" \
        "$OVERREP" \
        "$ADAPTER" >> "$OUT"

done

echo "QC summary created:"
echo "$OUT"
