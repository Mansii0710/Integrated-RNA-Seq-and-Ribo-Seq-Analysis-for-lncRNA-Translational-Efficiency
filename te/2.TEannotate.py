#!/usr/bin/env python3
import pandas as pd
import re

# --- Path to your GTF file ---
gtf_file = "/scratch/mallya/mansip/RNAseq/results/lnctpm_results/lncgenome/gencode_lncRNA.filtered.gtf"

# Function to parse attributes column
def parse_attributes(attr_string):
    attrs = {}
    for match in re.finditer(r'(\S+) "([^"]+)"', attr_string):
        attrs[match.group(1)] = match.group(2)
    return attrs

# Read GTF and extract annotation
records = []

with open(gtf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        
        if fields[2] != "transcript":
            continue
        
        attrs = parse_attributes(fields[8])
        
        records.append({
            "Transcript": attrs.get("transcript_id"),
            "gene_id": attrs.get("gene_id"),
            "gene_name": attrs.get("gene_name"),
            "gene_type": attrs.get("gene_type")
        })

gtf_df = pd.DataFrame(records).drop_duplicates()
print("Loaded transcripts from GTF:", gtf_df.shape)

# Load TE file
te_df = pd.read_csv("lncRNA_avg_TE_filtered.csv")

# Merge annotation
annotated = te_df.merge(gtf_df, on="Transcript", how="left")

# Save
annotated.to_csv("lncRNA_TE_annotated.csv", index=False)

print("Annotated TE saved to lncRNA_TE_annotated.csv")
print("Missing annotations:", annotated['gene_name'].isna().sum())
