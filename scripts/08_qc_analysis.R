# Sample-level QC for DE analysis (PCA + hierarchical clustering)
# 1. Load packages and the dds object from script 07
# pak::pkg_install(c("tidyverse", "DESeq2", "pheatmap", "RColorBrewer"))
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

# DESeqDataSet built in 07_count_normalization.R
dds <- read_rds("outputs/tables/dds.rds")

# 2. Regularized-log transform (for visualization only)
# blind = TRUE: transform without using the sample groups, so QC stays unbiased.
# (We have 16 samples, so rlog is fine; for >~20 samples prefer vst().)
rld <- rlog(dds, blind = TRUE)

# 3. Principal Component Analysis (PCA)
# Colour by our condition of interest (treatment). Replicates should group.
p_pca_cond <- plotPCA(rld, intgroup = "condition")
p_pca_cond

ggsave("outputs/figures/pca_condition.png", p_pca_cond,
       width = 7, height = 5, dpi = 300)

# Colour by donor (cell_line) to see the paired/batch structure we control for.
p_pca_line <- plotPCA(rld, intgroup = "cell_line")
p_pca_line
ggsave("outputs/figures/pca_cell_line.png", p_pca_line,
       width = 7, height = 5, dpi = 300)

# Save the PCA coordinates too
pca_data <- plotPCA(rld, intgroup = c("condition", "cell_line"),
                    returnData = TRUE)
write.csv(pca_data, "outputs/tables/qc_pca_data.csv", row.names = FALSE)

# 4. Hierarchical clustering heatmap (sample-sample correlation)
# Extract the rlog matrix, then correlate every pair of samples.
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)          # pairwise correlations (expect > 0.99)

# Metadata used to draw the coloured annotation bars on the heatmap
meta <- as.data.frame(colData(dds)[, c("condition", "cell_line")])

heat_colors <- brewer.pal(6, "Blues")
pheatmap(
  rld_cor,
  annotation    = meta,
  color         = heat_colors,
  border_color  = NA,
  fontsize      = 10,
  fontsize_row  = 10,
  filename      = "outputs/figures/sample_correlation_heatmap.png"
)