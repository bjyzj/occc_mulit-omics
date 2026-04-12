

library(tidyverse)
library(ggplot2)

protDE <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/proteomics_results/DE_OCCC_vs_ccRCC_proteomics_qn.csv")

sum(protDE$adj.P.Val < 0.05, na.rm = TRUE)
hist(protDE$logFC, breaks = 50)
hist(protDE$adj.P.Val, breaks = 50)

summary(protDE$logFC)

fc_thresh  <- 0.3 # 2-fold change threshold of +-1 in log2 scale
fdr_thresh <- 0.05 # FDR cutoff

volcano_df <- protDE %>% mutate(Significant = case_when(
  adj.P.Val < fdr_thresh & logFC > fc_thresh  ~ "Up",
  adj.P.Val < fdr_thresh & logFC < -fc_thresh ~ "Down", TRUE ~ "Not Significant"))
table(volcano_df$Significant)

ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "deeppink4",
                                "Down" = "dodgerblue4", 
                                "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed") +
  labs(x = "log2 Fold Change", y = "-log10(FDR)") +
  theme_classic()

#################


volcano_df <- protDE %>% mutate(Significant = case_when(
  adj.P.Val < fdr_thresh & logFC > fc_thresh  ~ "Up",
  adj.P.Val < fdr_thresh & logFC < -fc_thresh ~ "Down",
  TRUE ~ "Not Significant"
))

v3 <- ggplot(volcano_df, aes(x = logFC,
                             y = -log10(adj.P.Val),
                             color = Significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Up" = "deeppink4",
               "Down" = "dodgerblue4",
               "Not Significant" = "grey70"),
    labels = c(
      "Up" = "Upregulated\n(FDR < 0.05, log2FC > 0.3)",
      "Down" = "Downregulated\n(FDR < 0.05, log2FC < -0.3)",
      "Not Significant" = "Not significant"
    )
  ) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh),
             linetype = "dashed") +
  labs(
    x = "log2 Fold Change",
    y = expression(-log[10]("FDR"))
  ) +
  guides(color = guide_legend(title = NULL)) +
  theme_classic() +
  theme(
    legend.position = c(0.21, 0.80),
    legend.background = element_rect(
      fill = alpha("white", 0.7),
      color = "grey80",
      linewidth = 0.3
    ),
    legend.key.height = unit(0.8, "cm"),
    legend.text = element_text(size = 8)
  )

v3



