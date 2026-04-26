#!/usr/bin/env bash
#SBATCH --job-name=star_riboseq
#SBATCH --output=star_riboseq_%j.out
#SBATCH --error=star_riboseq_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --partition=standard

set -euo pipefail

# ==========================
# Activate conda environment
# ==========================
CONDA_BASE="$HOME/miniconda3"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate ribo_seq_env

echo "[INFO] Activated environment: $(conda info --envs | grep '*' )"
echo "[STEP] Aligning Ribo-seq reads with STAR..."

# ==========================
# PATHS
# ==========================
SCRATCH="/scratch/mallya/mansip"
FILTERED="$SCRATCH/ribseq2/results/rRNA_removal"
STAR_INDEX="/home/mansip/genome/star_index"
ALIGNMENTS="$SCRATCH/ribseq2/results/star_alignments"

mkdir -p "$ALIGNMENTS"

# ==========================
# MAIN STAR ALIGNMENT LOOP
# ==========================
for fq in "$FILTERED"/*_non_rRNA.fastq*; do

  sample_name=$(basename "$fq" | grep -oE 'SRR[0-9]+')
  echo "[INFO] Mapping $sample_name with STAR"

  STAR \
    --runThreadN $SLURM_CPUS_PER_TASK \
    --genomeDir "$STAR_INDEX" \
    --readFilesIn "$fq" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$ALIGNMENTS/${sample_name}_" \
    --outSAMtype BAM Unsorted \
    --quantMode TranscriptomeSAM \
    --outFilterMultimapNmax 1 \
    --alignEndsType EndToEnd \
    --alignIntronMax 1 \
    --outSAMattributes NH HI NM MD

  echo "[INFO] Completed STAR mapping for $sample_name"

  # ========= Sort genomic BAM =========
  samtools sort -@ $SLURM_CPUS_PER_TASK \
    -o "${ALIGNMENTS}/${sample_name}_Aligned.sortedByCoord.out.bam" \
    "${ALIGNMENTS}/${sample_name}_Aligned.out.bam"

  samtools index "${ALIGNMENTS}/${sample_name}_Aligned.sortedByCoord.out.bam"
  echo "[INFO] Sorted and indexed genome BAM for $sample_name"

  # ========= Sort transcriptome BAM =========
  samtools sort -@ $SLURM_CPUS_PER_TASK \
    -o "${ALIGNMENTS}/${sample_name}_Transcriptome.out.bam" \
    "${ALIGNMENTS}/${sample_name}_Aligned.toTranscriptome.out.bam"

  samtools index "${ALIGNMENTS}/${sample_name}_Transcriptome.out.bam"
  echo "[INFO] Sorted and indexed transcriptome BAM for $sample_name"

  # ========= Remove unsorted intermediate BAMs =========
  rm -f "${ALIGNMENTS}/${sample_name}_Aligned.out.bam"
  rm -f "${ALIGNMENTS}/${sample_name}_Aligned.toTranscriptome.out.bam"

done

echo "STAR alignment, sorting, renaming, and indexing complete."
