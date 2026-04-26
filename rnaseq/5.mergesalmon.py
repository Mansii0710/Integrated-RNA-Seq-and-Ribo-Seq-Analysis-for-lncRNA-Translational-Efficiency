#!/usr/bin/env python3
import os
import sys
import pandas as pd
from pathlib import Path

if len(sys.argv) < 3:
    print("Usage: python merge_tpm.py <salmon_results_dir> <output.tsv>")
    sys.exit(1)

results_dir = Path(sys.argv[1])
out_file = sys.argv[2]

dfs = []
sample_names = []

# Search for all quant.sf files
for q in sorted(results_dir.rglob("quant.sf")):
    sample = q.parent.name
    df = pd.read_csv(q, sep="\t", usecols=["Name", "TPM"])
    df = df.rename(columns={"TPM": sample})
    dfs.append(df.set_index("Name"))
    sample_names.append(sample)

# Merge all samples by transcript ID
merged = pd.concat(dfs, axis=1).fillna(0)

# Output
merged.to_csv(out_file, sep="\t")

print("Merged TPM matrix created:", out_file)
print("Shape:", merged.shape)

