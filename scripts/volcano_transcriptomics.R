
library(tidyverse)
library(ggrepel)
library(ggplot2)

DE_OCvsRC <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/transcriptomics_results/DE_results_OCCC_vs_ccRCC_ratio.csv")
DE_OCvsGTExOV <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/transcriptomics_results/DE_results_OCCC_vs_GTEx_Ovary_ratio.csv")

###########
# log2((OCCC / GTEx Ovary) / (ccRCC / GTEx Renal Cortex))
############
# fix zero adj p val to avoid inifinite log10 values
DE_OCvsRC$adj.P.Val[DE_OCvsRC$adj.P.Val == 0] <- 1e-300

# thresholds
fc_thresh  <- 1 # 2-fold change threshold of +-1 in log2 scale
fdr_thresh <- 0.05 # FDR cutoff

# significance
DE_OCvsRC$Significance <- "Not Significant"
DE_OCvsRC$Significance[DE_OCvsRC$adj.P.Val < fdr_thresh & DE_OCvsRC$logFC > fc_thresh] <- "Up"
DE_OCvsRC$Significance[DE_OCvsRC$adj.P.Val < fdr_thresh & DE_OCvsRC$logFC < -fc_thresh] <- "Down"



# Volcano plot
v1 <- ggplot(DE_OCvsRC, aes(x = logFC,
                            y = -log10(adj.P.Val),
                            color = Significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Up" = "red4",
                                "Down" = "navyblue",
                                "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh),
             linetype = "dashed") +
  labs(x = "log2 Fold Change", y = "-log10(FDR)") +
  theme_classic()

v1

#####################

v1 <- ggplot(DE_OCvsRC, aes(x = logFC,
                            y = -log10(adj.P.Val),
                            color = Significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("Up" = "red4",
               "Down" = "navyblue",
               "Not Significant" = "grey70"),
    labels = c(
      "Upregulated\n(FDR < 0.05, log2FC > 1)",
      "Downregulated\n(FDR < 0.05, log2FC < -1)",
      "Not significant"
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
    legend.position = c(0.82, 0.15),
    legend.background = element_rect(
      fill = alpha("white", 0.7),
      color = "grey80",
      linewidth = 0.3
    ),
    legend.key.height = unit(0.8, "cm"),
    legend.text = element_text(size = 8)
  )

v1
# the genes more upregulated/downregulated in OCCC relative to its normal tissues than in 
# ccRCC relative to its normal tissues -> specific or not specific to the cancer tissues


########
# OCCC vs GTEx
########

# fix zero adj p val to avoid inifinite log10 values
DE_OCvsGTExOV$adj.P.Val <- pmax(DE_OCvsGTExOV$adj.P.Val, 1e-300)

# thresholds
fc_thresh  <- 1 # 2-fold change threshold of +-1 in log2 scale
fdr_thresh <- 0.05 # FDR cutoff

# significance
DE_OCvsGTExOV$Significance <- "Not Significant"
DE_OCvsGTExOV$Significance[DE_OCvsGTExOV$adj.P.Val < fdr_thresh & DE_OCvsGTExOV$logFC > fc_thresh] <- "Up"
DE_OCvsGTExOV$Significance[DE_OCvsGTExOV$adj.P.Val < fdr_thresh & DE_OCvsGTExOV$logFC < -fc_thresh] <- "Down"



# Volcano plot
v2 <- ggplot(DE_OCvsGTExOV, aes(x = logFC,
                                y = -log10(adj.P.Val),
                                color = Significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Up" = "deeppink4",
                                "Down" = "dodgerblue4",
                                "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh),
             linetype = "dashed") +
  labs(x = "log2 Fold Change", y = "-log10(FDR)") +
  theme_classic()

v2


##################
v2 <- ggplot(DE_OCvsGTExOV, aes(x = logFC,
                            y = -log10(adj.P.Val),
                            color = Significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("Up" = "red4",
               "Down" = "navyblue",
               "Not Significant" = "grey70"),
    labels = c(
      "Upregulated\n(FDR < 0.05, log2FC > 1)",
      "Downregulated\n(FDR < 0.05, log2FC < -1)",
      "Not significant"
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
    legend.position = c(0.21, 0.15),
    legend.background = element_rect(
      fill = alpha("white", 0.7),
      color = "grey80",
      linewidth = 0.3
    ),
    legend.key.height = unit(0.8, "cm"),
    legend.text = element_text(size = 8)
  )

v2







