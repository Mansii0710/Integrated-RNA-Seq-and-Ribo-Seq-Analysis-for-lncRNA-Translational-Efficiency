#!/bin/bash
set -euo pipefail

SCRATCH="/scratch/mallya/mansip/ribseq2"  
RAW="$SCRATCH/rawfastq"
QC="$SCRATCH/results/qc/preqc"

mkdir -p "$QC"

echo "[STEP 3] Running FastQC on all FASTQ.gz files in $RAW"

# Find all fastq.gz files recursively in RAW folder
find "$RAW" -type f -name "*.fastq.gz" | while read -r filepath; do
  filename=$(basename "$filepath")
  echo "Processing file: $filename"
  fastqc "$filepath" -o "$QC" --quiet --threads 8
done

echo "[STEP 3] Running MultiQC on FastQC results"
multiqc "$QC" -o "$QC" --force --quiet

echo "[STEP 3] Done."
