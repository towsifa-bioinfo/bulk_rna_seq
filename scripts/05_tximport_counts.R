# RNA-seq Analysis in R (Part 1): tximport -> Gene-level Count Table
# 1. Install and load packages
# Run the install line ONCE, then comment it out.

pak::pkg_install(c("tidyverse", "tximport", "EnsDb.Hsapiens.v86","AnnotationDbi"))
#load packages
library(tidyverse)
library(tximport)
library(EnsDb.Hsapiens.v86)

# 2. Locate the Salmon quant files
# Each sample has its own folder under outputs/salmon_out/, containing quant.sf
samples <- list.dirs("outputs/salmon_out", recursive = FALSE, full.names = FALSE)
samples

# Build the full path to each sample's quant.sf and name it by sample
quant_files <- file.path("outputs/salmon_out", samples, "quant.sf")
names(quant_files) <- samples
print(quant_files)

# Safety check: every file must exist (all should print TRUE)
file.exists(quant_files) 

# 3. Build the transcript -> gene mapping (tx2gene)
# What are the columns in the database?
EnsDb.Hsapiens.v86
columns(EnsDb.Hsapiens.v86)
keys(EnsDb.Hsapiens.v86)

# Salmon counts transcripts; we want genes. 
# This table tells tximport which transcript belongs to which gene symbol.
# Get the mapping from transcript IDs to gene symbols 
tx2gene <- AnnotationDbi::select(
  EnsDb.Hsapiens.v86,
  keys    = keys(EnsDb.Hsapiens.v86),
  columns = c("TXID", "SYMBOL")
)

# Keep only transcript ID + symbol (drop the gene ID column)
tx2gene <- dplyr::select(tx2gene, TXID, SYMBOL)
head(tx2gene)

# 4. Import and summarize to gene level
txi <- tximport(
  files               = quant_files,
  type                = "salmon",
  tx2gene             = tx2gene,
  ignoreTxVersion     = TRUE,
  countsFromAbundance = "lengthScaledTPM"
)

class(txi)          # "list"
head(txi$counts)    # raw gene-level counts
head(txi$abundance) # gene-level TPM

# 5. Save the count table (the reusable product of this step)
# Raw counts (row.names = gene symbols) -> keep the gene column when writing CSV
raw_counts <- txi$counts
write.csv(raw_counts, "outputs/tables/raw_counts.csv", row.names = TRUE)
write_rds(raw_counts, "outputs/tables/raw_counts.rds")

# TPM (normalized abundance)
tpm_counts <- txi$abundance
write.csv(tpm_counts, "outputs/tables/tpm_counts.csv", row.names = TRUE)
write_rds(tpm_counts, "outputs/tables/tpm_counts.rds")

# Save the whole tximport object. The DESeq2 step needs it (not just the CSV),
# because DESeqDataSetFromTximport() uses the transcript-length information.
write_rds(txi, "outputs/tables/txi.rds")



