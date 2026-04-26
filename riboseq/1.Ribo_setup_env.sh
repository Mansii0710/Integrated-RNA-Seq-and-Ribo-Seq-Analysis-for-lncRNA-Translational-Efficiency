#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# Ribo-seq Environment Setup Script
# =============================================================================
# This script installs Miniconda (if needed), configures conda/mamba channels,
# and creates a reproducible ribo_seq_env with all tools required for
# ribosome profiling analysis (FASTQC, Cutadapt, TrimGalore, STAR, UMI-tools, etc.)
# Designed for HPC users; breakpoints allow troubleshooting.

# ----------------------
# 1. Variables
# ----------------------
INSTALL_DIR="$HOME/miniconda3"
ENV_NAME="ribo_seq_env"

# ----------------------
# 2. Install Miniconda
# ----------------------
if [[ ! -d "$INSTALL_DIR" ]]; then
    echo "[STEP 2] Installing Miniconda..."
    wget -O miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash miniconda.sh -b -p "$INSTALL_DIR"
    rm miniconda.sh
fi
source "$INSTALL_DIR/etc/profile.d/conda.sh"

# ----------------------
# 3. Configure Conda Channels
# ----------------------
echo "[STEP 3] Configuring conda channels..."
conda config --remove-key channels &>/dev/null || true
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict

# ----------------------
# 4. Install Mamba
# ----------------------
echo "[STEP 4] Installing mamba..."
conda install -y -n base -c conda-forge mamba

# =======================
# 5. Create Ribo-seq environment
# =======================
echo "[STEP 5] Creating conda environment: $ENV_NAME"
if ! conda env list | grep -qE "^$ENV_NAME[[:space:]]"; then
    mamba create -y -n $ENV_NAME \
    python=3.9 \
    fastqc=0.12.1 \
    multiqc=1.15 \
    cutadapt=4.5 \
    trim-galore=0.6.10 \
    bowtie2=2.5.1 \
    star=2.7.10b \
    samtools=1.17 \
    umi_tools=1.1.4 \
    salmon=1.10.1 \
    fq=0.11.0 \
    sortmerna=4.3.4 \
    -c conda-forge \
    -c bioconda \
    --channel-priority flexible

else
    echo "[STEP 5] Environment '$ENV_NAME' exists, skipping creation."
fi

# =======================
# 6. Activate environment
# =======================
echo "[STEP 6] Activating environment..."
conda activate $ENV_NAME
# STEP 7: Install quicksect via conda, then pip-install ribotricer without deps
echo "[STEP 7] Installing quicksect from conda-forge..."
mamba install -y -n $ENV_NAME quicksect -c conda-forge

echo "[STEP 8] Installing Ribo-TISH and Ribotricer via pip (no deps)..."
pip install --no-deps ribotish ribotricer

conda list --name $ENV_NAME | grep -E 'fastqc|cutadapt|trim-galore|bowtie2|star|samtools|umi_tools|sortmerna|salmon|ribotish|ribotricer'

echo "[STEP 9] Environment setup complete."
