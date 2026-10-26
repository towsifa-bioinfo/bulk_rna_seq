# Count Normalization with DESeq2 (median of ratios)
# 1. Load packages and inputs
# pak::pkg_install(c("tidyverse", "DESeq2"))
library(tidyverse)
library(DESeq2)

# tximport object from 02a_tximport_counts.R
txi <- read_rds("outputs/tables/txi.rds")
samples <- colnames(txi$counts)

# 2. Build metadata from the real design (SraRunTable), aligned to counts
run_table <- read.csv("inputs/SraRunTable.csv", check.names = FALSE)

meta <- run_table |> 
  dplyr::select(Run, treatment, cell_line)  |> 
  dplyr::distinct()

# CRITICAL: rows of metadata must be in the SAME ORDER as the count columns.
meta <- meta[match(samples, meta$Run), ]

col_data <- data.frame(
  row.names = meta$Run,
  condition = factor(meta$treatment),
  cell_line = factor(meta$cell_line)
)
col_data$condition <- relevel(col_data$condition, ref = "Untreated")

# DESeq2 errors out if names/order don't match - verify explicitly.
all(colnames(txi$counts) == rownames(col_data))
stopifnot(all(colnames(txi$counts) == rownames(col_data)))

dds <- DESeqDataSetFromTximport(
  txi     = txi,
  colData = col_data,
  design  = ~ cell_line + condition
)# 3. Create the DESeqDataSet
# Paired design: control for donor (cell_line), test the treatment effect.


head(counts(dds))   # raw count slot

# 4. Estimate size factors (median of ratios) and get normalized counts
dds <- estimateSizeFactors(dds)
dds

# Size factors should sit around 1; large swings hint at outlier samples.
sizeFactors(dds)
counts(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# 5. Save outputs for downstream visualization / QC
write.csv(normalized_counts, "outputs/tables/normalized_counts.csv", row.names = TRUE)
write_rds(normalized_counts, "outputs/tables/normalized_counts.rds")

# Save the dds object too, so later steps can reuse it.
write_rds(dds, "outputs/tables/dds.rds")