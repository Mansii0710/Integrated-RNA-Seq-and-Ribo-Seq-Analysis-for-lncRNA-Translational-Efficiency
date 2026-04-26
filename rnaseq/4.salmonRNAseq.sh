#!/bin/bash
#SBATCH --job-name=lncRNA_PE_full
#SBATCH --output=lncRNA_PE_full_%j.out
#SBATCH --error=lncRNA_PE_full_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --partition=standard

set -euo pipefail

##############################
# USER PATHS
##############################
GTF="/scratch/mallya/mansip/RNAseq/results/lnctpm_results/lncgenome/gencode.v48.long_noncoding_RNAs.gtf"
GENOME="/home/mansip/genome/GRCh38.primary_assembly.genome.fa"

FILTERED_GTF="/scratch/mallya/mansip/RNAseq/lncgtf_Salmon/lnc_filtered_PE.gtf"
FASTA="/scratch/mallya/mansip/RNAseq/lncgtf_Salmon/lncRNA.filtered_PE.fa"
INDEX="/scratch/mallya/mansip/RNAseq/lncgtf_Salmon/salmon_index_lncRNA_PE2"

TRIMMED="/scratch/mallya/mansip/RNAseq/trimmed"
OUTDIR="/scratch/mallya/mansip/lncgtf_Salmon/lncRNA_PE_quant"
mkdir -p $OUTDIR

##############################
# ENVIRONMENT
##############################

echo "=========================="
echo "STEP 1 FILTER LncRNA GTF"
echo "=========================="

# Keep primary chromosomes only
grep -E '^chr[0-9XYM]' $GTF > primary_lnc_PE.gtf

# Remove unwanted transcript types
awk '
/transcript_type "retained_intron"/ {next}
 /transcript_type "sense_intronic"/ {next}
 /transcript_type "antisense_intronic"/ {next}
 tolower($0) ~ /intron/ {next}
 {print}
' primary_lnc_PE.gtf > $FILTERED_GTF

echo "Filtered GTF saved as $FILTERED_GTF"


echo "=============================="
echo "STEP 2 Extract FASTA (gffread)"
echo "=============================="

gffread $FILTERED_GTF -g $GENOME -w $FASTA
echo "FASTA generated: $FASTA"


echo "=========================="
echo "STEP 3 Build Salmon Index"
echo "=========================="

rm -rf $INDEX   # Important fix
echo "Old index removed."

salmon index \
    -t $FASTA \
    -i $INDEX \
    --threads 16 \
    -k 31 \
    --gencode \

echo "Index created at:$INDEX"



echo "========================================"
echo "STEP 4 Salmon Quantification (Paired)"
echo "========================================"

cd $TRIMMED

for R1 in *_1_val_1.fq.gz
do
    R2="${R1/_1_val_1/_2_val_2}"
    SAMPLE=$(basename "$R1" "_1_val_1.fq.gz")

    echo "Quantifying sample: $SAMPLE"

    salmon quant \
        -i $INDEX \
        -l A \
        -1 "$R1" \
        -2 "$R2" \
        -p 16 \
        --validateMappings \
        --discardOrphansQuasi \
        -o "${OUTDIR}/${SAMPLE}"
done

echo "==========================="
echo "PIPELINE COMPLETE (PE MODE)"
echo "==========================="