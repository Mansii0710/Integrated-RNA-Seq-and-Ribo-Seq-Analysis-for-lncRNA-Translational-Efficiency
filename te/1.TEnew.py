#!/usr/bin/env python3
import pandas as pd

# Parameters
rna_file = "merge_RNA_TPM.csv"
ribo_file = "ribo_tpm_lncRNA_ONLY.csv"
min_rna_tpm = 1       # minimum average RNA TPM to keep transcript
min_ribo_tpm = 0.5    # minimum average Ribo TPM to keep transcript
pseudo = 1         # pseudo-count to avoid divide-by-zero

# Load RNA-seq TPMs
rna = pd.read_csv(rna_file, sep="\t", index_col=0)

# Load Ribo-seq TPMs
ribo = pd.read_csv(ribo_file, index_col=0)

# Drop extra columns from Ribo-seq if present
for col in ["tx_len", "biotype"]:
    if col in ribo.columns:
        ribo = ribo.drop(columns=[col])

# Keep only transcripts present in both
common = rna.index.intersection(ribo.index)
rna = rna.loc[common]
ribo = ribo.loc[common]

# Calculate average TPM per transcript
rna_avg = rna.mean(axis=1)
ribo_avg = ribo.mean(axis=1)

# Filter lowly expressed transcripts
keep = (rna_avg >= min_rna_tpm) & (ribo_avg >= min_ribo_tpm)
rna_avg = rna_avg[keep]
ribo_avg = ribo_avg[keep]

# Calculate TE per transcript
te = (ribo_avg + pseudo) / (rna_avg + pseudo)

# Convert to DataFrame and sort by TE
te_df = pd.DataFrame({
    "Transcript": te.index,
    "RNA_avg_TPM": rna_avg,
    "Ribo_avg_TPM": ribo_avg,
    "TE": te
}).sort_values("TE", ascending=False)

# Save to CSV
te_df.to_csv("lncRNA_avg_TE_filtered.csv", index=False)

print("Filtered TE per transcript saved to lncRNA_avg_TE_filtered.csv")
print("Number of transcripts after filtering:", te_df.shape[0])
