#!/bin/bash
#SBATCH --job-name=rnaseq_align_resume
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=24:00:00
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

# ===== 2?? Define paths =====

# Trimmed FASTQs directory (use directly!)
TRIM_HOME="/scratch/mallya/mansip/RNAseq/trimmed"

# Scratch folders for storing results
SCRATCH="/scratch/mallya/mansip/RNAseq/results"

ALIGN="$SCRATCH/star_align"
STRINGTIE="$SCRATCH/stringtie"
TRANSCRIPTS="$SCRATCH/transcripts"

# Genome / reference
GENOME_DIR="$HOME/genome"
REF_FASTA="$GENOME_DIR/GRCh38.primary_assembly.genome.fa"
REF_GTF="$GENOME_DIR/gencode.v48.annotation.gtf"
STAR_INDEX="$GENOME_DIR/star_index"

THREADS=$SLURM_CPUS_PER_TASK

mkdir -p "$ALIGN" "$STRINGTIE" "$TRANSCRIPTS"

echo "[INFO] Using trimmed FASTQs directly from: $TRIM_HOME"
echo "[INFO] Storing alignment & transcript results in: $SCRATCH"

# ===== 3?? No copying step � Directly process files =====

cd "$TRIM_HOME"

# ===== 4?? STAR + StringTie =====
for R1 in *_1_val_1.fq.gz; do
    SAMPLE=$(basename "$R1" _1_val_1.fq.gz)
    R2="${SAMPLE}_2_val_2.fq.gz"
    BAM="$ALIGN/${SAMPLE}.bam"

    echo "[INFO] Processing sample: $SAMPLE"

    STAR --runThreadN $THREADS \
         --genomeDir "$STAR_INDEX" \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$ALIGN/${SAMPLE}_" \
         --outSAMtype BAM SortedByCoordinate

    mv "$ALIGN/${SAMPLE}_Aligned.sortedByCoord.out.bam" "$BAM"
    samtools index "$BAM"

    # StringTie Assembly
    GTF="$STRINGTIE/${SAMPLE}.gtf"
    stringtie "$BAM" -G "$REF_GTF" -p $THREADS -o "$GTF"

    # Known vs novel transcripts split
    awk '/transcript/ && /reference_id/ {keep=1} /transcript/ && !/reference_id/ {keep=0} keep || /reference_id/' \
        "$GTF" > "$STRINGTIE/${SAMPLE}_known.gtf"

    awk '/transcript/ && !/reference_id/ {keep=1} /transcript/ && /reference_id/ {keep=0} keep' \
        "$GTF" > "$STRINGTIE/${SAMPLE}_novel.gtf"

    # Extract transcript FASTA
    gffread "$GTF" -g "$REF_FASTA" -w "$TRANSCRIPTS/${SAMPLE}_all.fa"
    gffread "$STRINGTIE/${SAMPLE}_known.gtf" -g "$REF_FASTA" -w "$TRANSCRIPTS/${SAMPLE}_known.fa"
    gffread "$STRINGTIE/${SAMPLE}_novel.gtf" -g "$REF_FASTA" -w "$TRANSCRIPTS/${SAMPLE}_novel.fa"
done

echo "[INFO] Results stored in: $SCRATCH"
echo "[INFO] Alignment + Assembly completed successfully!"
