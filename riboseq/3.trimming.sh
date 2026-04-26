#!/bin/bash
#SBATCH --job-name=ribo_cutadapt_se
#SBATCH --output=ribo_cutadapt_se_%j.out
#SBATCH --error=ribo_cutadapt_se_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --partition=hm

# ======== SETUP CONDA ENVIRONMENT ========
CONDA_BASE="$HOME/miniconda3"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate ribo_seq_env

# ======== PATHS ========
RAW_DIR=/scratch/mallya/mansip/ribseq2/rawfastq
TRIMMED=/scratch/mallya/mansip/ribseq2/trimmed
QC=/scratch/mallya/mansip/ribseq2/results/qc/postqc

mkdir -p ${TRIMMED} ${QC}

# ======== ADAPTER ========
ADAPTOR=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

# ======== LOOP OVER ALL FASTQ FILES ========
for IN in ${RAW_DIR}/*.fastq.gz; do
    BASENAME=$(basename ${IN} .fastq.gz)

    # ======== RUN CUTADAPT ========
    cutadapt \
      -a ${ADAPTOR} \
      --trim-n -q 20 \
      --minimum-length 18 --maximum-length 32 \
      --cores ${SLURM_CPUS_PER_TASK} \
      -o ${TRIMMED}/${BASENAME}.trim.fastq.gz \
      ${IN}

    # ======== RUN FASTQC ========
    fastqc -t 4 -o ${QC} ${TRIMMED}/${BASENAME}.trim.fastq.gz

    echo "Finished trimming and QC for ${BASENAME} at $(date)"
done
