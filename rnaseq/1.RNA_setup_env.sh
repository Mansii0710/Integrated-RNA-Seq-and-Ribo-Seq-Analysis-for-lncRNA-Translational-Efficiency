#!/bin/bash
set -euo pipefail

ENV_NAME="rnaseq_env"

echo "[INFO] Setting up environment: $ENV_NAME"

# Add conda to PATH and initialize
export PATH="/home/mansip/miniconda3/bin:$PATH"
source "/home/mansip/miniconda3/etc/profile.d/conda.sh"

# Configure channels
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority flexible

# Install mamba if missing
if ! command -v mamba &>/dev/null; then
    echo "[INFO] Installing mamba into base environment..."
    conda install -y -n base -c conda-forge mamba
fi

# Create or update environment
if conda env list | grep -q "^$ENV_NAME "; then
    echo "[INFO] Environment $ENV_NAME already exists, updating..."
    mamba install -y -n $ENV_NAME \
        python=3.9 \
        star=2.7.10b \
        stringtie \
        fastqc \
        gffread \
        gffcompare \
        seqkit \
        samtools \
        htslib \
        wget \
        rsync \
        pigz \
        isa-l \
        dos2unix \
        trim-galore \
        cutadapt=5.1
else
    echo "[INFO] Creating new environment $ENV_NAME..."
    mamba create -y -n $ENV_NAME \
        python=3.9 \
        star=2.7.10b \
        stringtie \
        fastqc \
        gffread \
        gffcompare \
        seqkit \
        samtools \
        htslib \
        wget \
        rsync \
        pigz \
        isa-l \
        dos2unix \
        trim-galore \
        cutadapt=5.1
fi

# Activate environment
conda activate $ENV_NAME

# Fix library path
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}"

# Verify installations
echo "[INFO] Environment setup completed ✅"
echo "[INFO] Installed tool versions:"
STAR --version
stringtie --version
fastqc --version
gffread --version
gffcompare --version
seqkit version
samtools --version | head -n1
trim_galore --version
cutadapt --version
