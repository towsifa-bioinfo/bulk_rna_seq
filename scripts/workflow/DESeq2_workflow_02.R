# install package
pak::pkg_install(c("tidyverse", "recount", "airway", "DESeq2", "EnhancedVolcano"))

# load the packages
library(tidyverse)
library(recount)
library(DESeq2)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Problem: Effects of EWS-FLI1 KD in Ewing sarcoma
# RQ1: What is the effect of knocking down EWS-FLI1 in Ewing sarcoma?
# RQ2: How does it change cellular behavior?

# The dataset we use here is downloaded from recount2
# Find a project of interest (note: it's hit or miss, the database is not limitless)
project_info <- abstract_search(query = "Ewing sarcoma")

# select study
download_study("SRP015989")

# load the data
load("SRP015989/rse_gene.Rdata")

# examine the read counts
counts_data <- assay(rse_gene)

# examine the coldata
col_data <- as.data.frame(colData(rse_gene))

# fix the colData to give a column with the appropriate groups
rse_gene$condition <- c(rep("shCTR", 3), rep("ShEF1", 4))

# check colnames of count_data and rownames of col_data
colnames(counts_data)
rownames(col_data)
all(colnames(counts_data) %in% rownames(col_data))

# are they in same order
all(colnames(counts_data) == rownames(col_data))

# perform DESeq2 analysis using DESeqDataSet
dds <- DESeqDataSet(rse_gene, design = ~condition)

# normalize the data
rld <- rlog(dds)
plotPCA(rld)

# perform DEG analysis
dds <- DESeq(dds)

# save as results
res <- results(dds)
res

# exploring the results
summary(res)


# add gene annotation to the results
# 1. check gene ID
rownames(res)

ens_ids <- sub("\\..*", "", rownames(res))
res$symbol <- mapIds(org.Hs.eg.db, keys = ens_ids,
                     keytype = "ENSEMBL", column = "SYMBOL",
                     multiVals = "first")

res$entrez <- mapIds(org.Hs.eg.db, keys = ens_ids,
                     keytype = "ENSEMBL", column = "ENTREZID",
                     multiVals = "first")

res$genename <- mapIds(org.Hs.eg.db, keys = ens_ids,
                       keytype = "ENSEMBL", column = "GENENAME",
                       multiVals = "first")

# dds df
anno_results <- as.data.frame(res)

# Volcano Plot
# Default cutoffs are log2FC > |2| and adjusted P-value < 0.05
EnhancedVolcano(
  anno_results,
  x = "log2FoldChange",
  y = "padj",
  lab = anno_results$symbol,
  pCutoff = 0.01,
  FCcutoff = 2,
  pointSize = 1.8,
  labSize = 3.0,
  selectLab = c("SOX2", "CLSTN2", "CD44", "MYC", "ERG"),  
  drawConnectors = TRUE,
  title = "EWS-FLI1 Knockdown"
)