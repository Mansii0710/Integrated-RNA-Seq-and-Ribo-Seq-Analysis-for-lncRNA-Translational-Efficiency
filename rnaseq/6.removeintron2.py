import pandas as pd

# Input and output file paths
input_csv = "ribo_tpm_lncRNA_all.csv"
output_csv = "ribo_tpm_lncRNA_ONLY.csv"

# Load file
df = pd.read_csv(input_csv)

# Keep only rows where biotype == "lncRNA"
df_filtered = df[df["biotype"] == "lncRNA"]

# Save
df_filtered.to_csv(output_csv, index=False)

print(f"Done! Saved filtered file with only lncRNA biotype to: {output_csv}")
print(f"Original rows: {len(df)}, Filtered rows (lncRNA only): {len(df_filtered)}")