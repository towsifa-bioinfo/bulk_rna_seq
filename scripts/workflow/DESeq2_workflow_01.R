# install package
pak::pkg_install(c("tidyverse","airway","DESeq2","EnhancedVolcano","org.Hs.eg.db","AnnotationDbi"))

#package for gene anotations
library(org.Hs.eg.db)
library(AnnotationDbi)


#load the packages 
library(tidyverse)
library(airway)
library(DESeq2)
library(EnhancedVolcano)
library(tidyplots)

# get the info about data
?airway
data("airway")

# counts table
counts_data <- assay(airway)

# metadata (col_data)
col_data <- as.data.frame(colData(airway))

#making sure the row names in `col_data` (metadata) matches to columns
rownames(col_data)
colnames(counts_data)
all(rownames(col_data)%in% colnames(counts_data))

#are they in same order
all(rownames(col_data) == colnames(counts_data))

#construct a `DESeqSetFromMatrix` data object
dds <- DESeqDataSetFromMatrix( 
  countData = counts_data,
  colData = col_data,
  # condition
  design = ~dex
)
dds 
# pre-filtering removing rows with low gene counts
# keep rows that have at least 10 reads total
rowSums(counts(dds))
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep, ]

# refernce category ~ set the factor level
dds$dex
class(dds$dex)
dds$dex <- relevel(dds$dex, ref = "untrt")

# perform different gene expression(DEG)
dds <- DESeq(dds)

# check differential gene expression analysis results
res <- results(dds)
res
names(res)

#exploring results
summary(res)

#signficance level (alpha)
#working with alpha
#if alpha = 0.1 or 10% ~ 90% CI (default)
summary(res)

# if alpha = 0.01 or 1% ~ 99% CI 
res1 <- results(dds, alpha = 0.01)
summary(res1)

#if alpha = 0.05 or 5% ~ 95% CI (most used)
res2 <- results(dds, alpha = 0.05)
summary(res2)

#results name
resultsNames(dds)

# contrast
contrast_res <- results(dds, contrast  = c("dex", "trt","untrt") )
contrast_res

# dds results into df
res_df <- as.data.frame(res2)

# add gene annotation to the results
# 1.check gene ID
rownames(res2)
# 2.column names of db
AnnotationDbi::columns(org.Hs.eg.db)

anno_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  rownames(res2),
  columns = c("ENSEMBL","ENTREZID", "SYMBOL","GENENAME" ),
  keytype = "ENSEMBL"
)
#check names
names(anno_df)

# dds df
res_df <- as.data.frame(res2)
res_df <- cbind(ENSEMBL = rownames(res_df), res_df)
rownames(res_df) <- NULL

# join annotation table to results
anno_results <- left_join(res_df, anno_df, by = "ENSEMBL")

#volcano plot
#Default cutoffs are log2FC > |2| andadjusted P-value < 0.05
EnhancedVolcano(
  anno_results,
  x = "log2FoldChange",
  y = "padj",
  lab = anno_results$SYMBOL,
  pCutoff = 0.01,
  FCcutoff = 2,
  pointSize = 1.5,
  labSize = 3.0,
  title = "Untreated vs. Treated"
)

