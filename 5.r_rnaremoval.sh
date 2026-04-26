#!/bin/bash
#SBATCH --job-name=sortmerna_human_rrna_removal
#SBATCH --output=sortmerna_%j.out
#SBATCH --error=sortmerna_%j.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --partition=standard

# ======== ACTIVATE CONDA SORTMERNA ENV ========
CONDA_BASE="$HOME/miniconda3"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate ribo_seq_env

# ======== DIRECTORIES ========
HOME_FASTA="/scratch/mallya/mansip/ribseq2/R_rna"          # rRNA reference fasta files
TRIMMED_DIR="/scratch/mallya/mansip/ribseq2/trimmed"       # cutadapt trimmed FASTQ
FILTERED_DIR="/scratch/mallya/mansip/ribseq2/results/rRNA_removal"
SCRATCH_RUN="/scratch/mallya/mansip/ribseq2/sortmerna_work"

mkdir -p "$FILTERED_DIR" "$SCRATCH_RUN"

# ======== rRNA Reference Database Files ========
RRNA_DB_5_8S="$HOME_FASTA/rfam-5.8s.fasta"
RRNA_DB_5S="$HOME_FASTA/rfam-5s.fasta"
RRNA_DB_18S="$HOME_FASTA/silva-euk-18s-human.fasta"
RRNA_DB_28S="$HOME_FASTA/silva-euk-28s-human.fasta"

# ======== GET ALL TRIMMED FASTQ FILES ========
samples=(${TRIMMED_DIR}/*.trim.fastq.gz)

echo "Found ${#samples[@]} trimmed FASTQ files."

for sample_path in "${samples[@]}"; do
    sample=$(basename "$sample_path")
    sample_prefix="${sample%.trim.fastq.gz}"

    echo "Processing sample: $sample"

    # Create unique SortMeRNA workdir
    WORKDIR="$SCRATCH_RUN/${sample_prefix}_work"
    mkdir -p "$WORKDIR"

    # ======== RUN SORTMERNA ========
    sortmerna \
        --ref "$RRNA_DB_5_8S" \
        --ref "$RRNA_DB_5S" \
        --ref "$RRNA_DB_18S" \
        --ref "$RRNA_DB_28S" \
        --reads "$sample_path" \
        --threads $SLURM_CPUS_PER_TASK \
        --workdir "$WORKDIR" \
        --aligned "$FILTERED_DIR/${sample_prefix}_rRNA.fastq" \
        --other "$FILTERED_DIR/${sample_prefix}_non_rRNA.fastq" \
        --fastx \
        --num_alignments 1 \
        -v

    echo "Finished: $sample"
done

echo " Human rRNA removal completed for ALL samples."
