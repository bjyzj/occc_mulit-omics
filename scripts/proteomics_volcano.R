

library(tidyverse)
library(ggplot2)

protDE <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/proteomics_results/DE_OCCC_vs_ccRCC_proteomics.csv")


fc_thresh  <- 1 # 2-fold change threshold of +-1 in log2 scale
fdr_thresh <- 0.05 # FDR cutoff

volcano_df <- protDE %>% mutate(Significant = case_when(
  adj.P.Val < fdr_thresh & logFC > fc_thresh  ~ "Up",
  adj.P.Val < fdr_thresh & logFC < -fc_thresh ~ "Down", TRUE ~ "Not Significant"))


ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "deeppink4",
                                "Down" = "dodgerblue4",
                                "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed") +
  labs(x = "log2 Fold Change", y = "-log10(FDR)") +
  theme_classic()




