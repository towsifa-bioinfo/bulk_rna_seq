# RNA-seq Count Data Exploration (concepts / QC before DESeq2)
# 1. Load packages and the count table from script 05
# pak::pkg_install("tidyverse")
library(tidyverse)

# Gene-level counts produced by 02a_tximport_counts.R
raw_counts <- readRDS("outputs/tables/raw_counts.rds")

# Round to whole numbers just for viewing / plotting
data <- as.data.frame(round(raw_counts))

# 2. Distribution of counts for a single sample
# Using SRR1039508 (an Untreated sample). RNA-seq counts are NOT normal:
# many low-count genes, a long right tail, and a huge dynamic range.
p_hist <- ggplot(raw_counts) +
  geom_histogram(aes(x = SRR1039508), bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
p_hist

ggsave("outputs/figures/count_distribution_SRR1039508.png",
       p_hist, width = 7, height = 5, dpi = 300)

# 3. Mean vs variance across the Untreated replicates
# The four Untreated replicates in our design
untreated <- c("SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520")

# For each gene, mean and variance across those four samples
mean_counts     <- apply(raw_counts[, untreated], 1, mean)
variance_counts <- apply(raw_counts[, untreated], 1, var)
df <- data.frame(mean_counts, variance_counts)

# Red line = where mean would equal variance (the Poisson assumption).
# Points sitting ABOVE the line show variance > mean (overdispersion),
# which is why we use the Negative Binomial (DESeq2), not Poisson.
p_mv <- ggplot(df) +
  geom_point(aes(x = mean_counts, y = variance_counts), alpha = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  xlab("Mean count (Untreated replicates)") +
  ylab("Variance (Untreated replicates)")

ggsave("outputs/figures/mean_vs_variance_untreated.png",
       p_mv, width = 7, height = 5, dpi = 300)