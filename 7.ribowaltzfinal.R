########################################
# Ribowaltz + lncRNA TPM Pipeline (All lncRNAs included) + QC
########################################

library(riboWaltz)
library(data.table)
library(dplyr)
library(ggplot2)
library(rtracklayer)

#--------------------------
# Paths
#--------------------------
bamfolder <- "/scratch/mallya/mansip/ribseq2/results/star_alignments/Transcriptomebam"
full_gtf  <- "/home/mansip/genome/gencode.v48.annotation.gtf"       
lnc_gtf   <- "/scratch/mallya/mansip/RNAseq/results/lnctpm_results/lncgenome/gencode.v48.long_noncoding_RNAs.gtf"
output_dir <- "/scratch/mallya/mansip/ribseq2/results/star_alignments/Transcriptomebam"
dir.create(output_dir, showWarnings = FALSE)

#--------------------------
# Load full transcriptome annotation (for P-site mapping)
#--------------------------
annotation <- create_annotation(gtfpath = full_gtf)
cat("Total transcripts in full GTF:", length(annotation$transcript), "\n")

#--------------------------
# List ONLY transcriptome BAM files
#--------------------------
bam_files <- list.files(
    bamfolder,
    pattern = "Transcriptome\\.out\\.bam$",
    full.names = TRUE
)

if(length(bam_files) == 0)
    stop("No Transcriptome.out.bam files found! Check path/pattern.")

cat("Using BAM files:\n")
print(bam_files)

#--------------------------
# Load reads
#--------------------------
reads_list <- bamtolist(
  bamfolder = bamfolder,
  annotation = annotation,
  transcript_align = TRUE,
  indel_threshold = 5,
  output_class = "list"
)

#--------------------------
# Compute P-site offsets
#--------------------------
psite_res <- psite(reads_list, plot = TRUE)
write.csv(psite_res, file.path(output_dir, "psite_offsets_all_samples.csv"), row.names = FALSE)

#--------------------------
# Flattened P-sites
#--------------------------
p_sites <- psite_info(reads_list, psite_res)
samples <- names(p_sites)
for (s in samples) {
  fwrite(p_sites[[s]], file.path(output_dir, paste0("psites_", s, ".csv")))
}

#--------------------------
# Aggregate transcript-level counts per sample
#--------------------------
ribo_counts <- lapply(names(p_sites), function(sample_name){
  dt <- p_sites[[sample_name]]
  if(nrow(dt) == 0) return(NULL)
  counts <- dt[, .(Ribo_count = .N), by = transcript]
  setnames(counts, "Ribo_count", sample_name)
  counts
})
names(ribo_counts) <- names(p_sites)
ribo_counts <- ribo_counts[!sapply(ribo_counts, is.null)]
if(length(ribo_counts) == 0) stop("No P-sites detected in any sample!")

# Merge all samples
ribo_matrix <- Reduce(function(x, y) merge(x, y, by="transcript", all=TRUE, sort=FALSE), ribo_counts)
ribo_matrix[is.na(ribo_matrix)] <- 0

#--------------------------
# Load lncRNA GTF for transcript lengths
#--------------------------
lnc_gtf_obj <- import(lnc_gtf)
lnc_tx_info <- lnc_gtf_obj[lnc_gtf_obj$type == "transcript"]
tx_lengths <- unique(data.table(
  tx_name = lnc_tx_info$transcript_id,
  tx_len  = width(lnc_tx_info),
  biotype = lnc_tx_info$transcript_type
), by = "tx_name")

#--------------------------
# Ensure all lncRNAs are included (even if no reads)
#--------------------------
all_lnc <- data.table(transcript = tx_lengths$tx_name)
ribo_matrix_lnc <- merge(all_lnc, ribo_matrix, by = "transcript", all.x = TRUE)
ribo_matrix_lnc[is.na(ribo_matrix_lnc)] <- 0

# Merge transcript lengths and biotype
ribo_matrix_lnc <- merge(
  ribo_matrix_lnc,
  tx_lengths,
  by.x = "transcript",
  by.y = "tx_name",
  all.x = TRUE,
  sort = FALSE
)

#--------------------------
# Convert counts to numeric
#--------------------------
sample_cols <- setdiff(colnames(ribo_matrix_lnc), c("tx_len","biotype","transcript"))
ribo_matrix_lnc[, (sample_cols) := lapply(.SD, as.numeric), .SDcols = sample_cols]

#--------------------------
# Compute TPM
#--------------------------
ribo_tpm <- ribo_matrix_lnc
for(col in sample_cols){
  rpk <- ribo_matrix_lnc[[col]] / (ribo_matrix_lnc$tx_len / 1000)
  if(sum(rpk) == 0){
    ribo_tpm[[col]] <- 0
  } else {
    ribo_tpm[[col]] <- (rpk / sum(rpk)) * 1e6
  }
}

#--------------------------
# Save TPM and counts
#--------------------------
fwrite(ribo_matrix_lnc, file.path(output_dir, "ribo_counts_lncRNA_all.csv"))
fwrite(ribo_tpm, file.path(output_dir, "ribo_tpm_lncRNA_all.csv"))
cat("\nPipeline completed! Counts and TPMs for all lncRNAs saved in:", output_dir, "\n")