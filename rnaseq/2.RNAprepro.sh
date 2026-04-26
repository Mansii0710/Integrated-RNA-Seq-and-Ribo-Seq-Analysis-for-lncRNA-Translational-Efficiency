#!/bin/bash
#SBATCH --job-name=rnaseq_preqc
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

set -euo pipefail

# ===== 1?? Load conda environment =====
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate rnaseq_env
    export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"
else
    echo "[ERROR] Conda environment not found!"
    exit 1
fi

# ===== 2?? Define paths (NO raw_fastq) =====
SCRATCH="/scratch/mallya/mansip/RNAseq"
RAW_HOME="$SCRATCH/rawfastq"          # Input FASTQs

# QC folder structure
QC_FASTQC_RAW="$SCRATCH/qc/fastqc/raw"
QC_FASTQC_TRIM="$SCRATCH/qc/fastqc/trimmed"
QC_MQC_RAW="$SCRATCH/qc/multiqc/raw"
QC_MQC_TRIM="$SCRATCH/qc/multiqc/trimmed"

TRIM="$SCRATCH/trimmed"

THREADS=$SLURM_CPUS_PER_TASK

# Create directories
mkdir -p "$QC_FASTQC_RAW" "$QC_FASTQC_TRIM" "$QC_MQC_RAW" "$QC_MQC_TRIM" "$TRIM"

echo "[INFO] Raw FASTQ directory: $RAW_HOME"

cd "$RAW_HOME"

# ===== 3?? FastQC on raw FASTQ =====
echo "[INFO] Running FastQC on raw reads..."
fastqc -t $THREADS -o "$QC_FASTQC_RAW" *.fastq.gz

# ===== 4?? TrimGalore =====
echo "[INFO] Running TrimGalore..."
for R1 in *_1.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="${SAMPLE}_2.fastq.gz"

    if [ ! -f "$R2" ]; then
        echo "[WARNING] Missing pair for $R1 � skipping!"
        continue
    fi

    trim_galore --paired --cores $THREADS -o "$TRIM" "$R1" "$R2"
done

# ===== 5?? FastQC on trimmed reads =====
echo "[INFO] FastQC on trimmed reads..."
fastqc -t $THREADS -o "$QC_FASTQC_TRIM" "$TRIM"/*_val_*.fq.gz

# ===== 6?? MultiQC =====
echo "[INFO] Running MultiQC..."
multiqc "$QC_FASTQC_RAW"   -o "$QC_MQC_RAW"   -n "multiqc_raw.html"
multiqc "$QC_FASTQC_TRIM"  -o "$QC_MQC_TRIM"  -n "multiqc_trimmed.html"

echo "[INFO] All processing complete! Output in: $SCRATCH"
