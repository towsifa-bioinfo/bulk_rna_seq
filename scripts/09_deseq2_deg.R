# RNA-seq Analysis in R (Part 2): DESeq2 Differential Expression
# 1. Install and load packages
pak::pkg_install(c("tidyverse", "DESeq2"))

library(tidyverse)
library(DESeq2)

# 2. Load the tximport object from step 04
txi <- read_rds("outputs/tables/txi.rds")
samples <- colnames(txi$counts)   # SRR ids, in the exact order DESeq2 will use
samples

# 3. Build metadata from the real design table
# SraRunTable.csv maps each SRR run to its treatment and cell line.
# This REPLACES the old hardcoded rep() metadata, which mislabeled samples.
run_table <- read.csv("inputs/SraRunTable.csv", check.names = FALSE)

# Keep the columns we need and key by the SRR run id
meta <- run_table  |> 
  dplyr::select(Run, treatment, cell_line)  |> 
  dplyr::distinct()

# CRITICAL: reorder metadata to match the column order of the counts.
# DESeq2 assumes row i of colData describes column i of the counts.
meta <- meta[match(samples, meta$Run), ]
identical(meta$Run, samples)
stopifnot(identical(meta$Run, samples))   # must line up exactly

col_data <- data.frame(
  row.names = meta$Run,
  condition = factor(meta$treatment),
  cell_line = factor(meta$cell_line)
)

# Make "Untreated" the reference level so results read as "X vs Untreated".
col_data$condition <- relevel(col_data$condition, ref = "Untreated")

col_data
table(col_data$condition)   # sanity check: 4 per treatment

# 4. Build the DESeq2 dataset
# The study is PAIRED by cell line (same 4 donors across treatments), so we
# control for cell_line and test the treatment effect: ~ cell_line + condition
dds <- DESeqDataSetFromTximport(
  txi     = txi,
  colData = col_data,
  design  = ~ cell_line + condition
)

# 5. Quality look: PCA
rlog_dds <- rlog(dds)
plotPCA(rlog_dds, intgroup = "condition")

pca_data <- plotPCA(rlog_dds, intgroup = c("condition", "cell_line"),
                    returnData = TRUE)
write.csv(pca_data, "outputs/tables/pca_data.csv", row.names = FALSE)
write_rds(pca_data, "outputs/tables/pca_data.rds")

# 6. Differential expression
dds <- DESeq(dds)
resultsNames(dds)   # shows which contrasts are available

# Main comparison for this study: Dexamethasone vs Untreated
resdf <- results(dds, contrast = c("condition", "Dexamethasone", "Untreated"))

plotMA(resdf)

# Tidy into a data frame with a gene column, then save

resdf <- as.data.frame(resdf)
resdf$gene <- rownames(resdf)
rownames(resdf) <- NULL
resdf <- dplyr::relocate(resdf, gene)

write.csv(resdf, "outputs/tables/DESeq2_results.csv", row.names = FALSE)
write_rds(resdf, "outputs/tables/DESeq2_results.rds")