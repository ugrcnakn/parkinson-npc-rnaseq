#!/bin/bash

# ==========================================
# Parkinson NPC RNA-seq
# Download + Validate + FASTQ + FastQC
# ==========================================

PROJECT="$HOME/bioinformatics/parkinson_npc_rnaseq"

ACCESSIONS="$PROJECT/data/metadata/remaining_accessions.txt"

RAW_DIR="$PROJECT/data/raw"
FASTQ_DIR="$PROJECT/data/raw/fastq"

QC_DIR="$PROJECT/results/fastqc/all_samples"
LOG_DIR="$PROJECT/logs/download_qc"

TMP_DIR="$PROJECT/tmp/fasterq"


# ------------------------------------------
# Create required directories
# ------------------------------------------

mkdir -p "$RAW_DIR"
mkdir -p "$FASTQ_DIR"
mkdir -p "$QC_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$TMP_DIR"


# ------------------------------------------
# Process each SRA accession
# ------------------------------------------

while IFS= read -r SRR || [[ -n "$SRR" ]]
do

    # Skip empty lines
    [[ -z "$SRR" ]] && continue

    LOG="$LOG_DIR/${SRR}.log"

    echo "==========================================" | tee "$LOG"
    echo "Starting: $SRR" | tee -a "$LOG"
    echo "Time: $(date)" | tee -a "$LOG"
    echo "==========================================" | tee -a "$LOG"


    # --------------------------------------
    # STEP 1 - Download SRA
    # --------------------------------------

    if [[ -s "$RAW_DIR/$SRR/$SRR.sra" ]]
    then
        echo "SRA already exists: $SRR" | tee -a "$LOG"
    else
        echo "Downloading SRA: $SRR" | tee -a "$LOG"

        if ! (
            cd "$RAW_DIR" &&
            prefetch "$SRR"
        ) >> "$LOG" 2>&1
        then
            echo "ERROR: prefetch failed for $SRR" | tee -a "$LOG"
            continue
        fi
    fi


    # --------------------------------------
    # STEP 2 - Validate SRA
    # --------------------------------------

    echo "Validating SRA: $SRR" | tee -a "$LOG"

    if ! vdb-validate "$RAW_DIR/$SRR" >> "$LOG" 2>&1
    then
        echo "ERROR: validation failed for $SRR" | tee -a "$LOG"
        continue
    fi


    # --------------------------------------
    # STEP 3 - Convert SRA to FASTQ
    # --------------------------------------

    if [[ -s "$FASTQ_DIR/${SRR}.fastq" ]]
    then
        echo "FASTQ already exists: $SRR" | tee -a "$LOG"
    else
        echo "Converting to FASTQ: $SRR" | tee -a "$LOG"

        mkdir -p "$TMP_DIR/$SRR"

        if ! fasterq-dump "$RAW_DIR/$SRR" \
            -O "$FASTQ_DIR" \
            -t "$TMP_DIR/$SRR" \
            -p >> "$LOG" 2>&1
        then
            echo "ERROR: fasterq-dump failed for $SRR" | tee -a "$LOG"
            continue
        fi
    fi


    # --------------------------------------
    # STEP 4 - Check FASTQ
    # --------------------------------------

    if [[ ! -s "$FASTQ_DIR/${SRR}.fastq" ]]
    then
        echo "ERROR: FASTQ file missing for $SRR" | tee -a "$LOG"
        continue
    fi


    # --------------------------------------
    # STEP 5 - Run FastQC
    # --------------------------------------

    if [[ -s "$QC_DIR/${SRR}_fastqc.html" ]]
    then
        echo "FastQC already exists: $SRR" | tee -a "$LOG"
    else
        echo "Running FastQC: $SRR" | tee -a "$LOG"

        if ! fastqc "$FASTQ_DIR/${SRR}.fastq" \
            -o "$QC_DIR" >> "$LOG" 2>&1
        then
            echo "ERROR: FastQC failed for $SRR" | tee -a "$LOG"
            continue
        fi
    fi


    # --------------------------------------
    # SAMPLE COMPLETE
    # --------------------------------------

    echo "Completed successfully: $SRR" | tee -a "$LOG"
    echo "Time: $(date)" | tee -a "$LOG"
    echo "" | tee -a "$LOG"

done < "$ACCESSIONS"

echo "All accessions processed."
