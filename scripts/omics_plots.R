
library(tidyverse)
library(ggplot2)
library(limma)
#library(ggVennDiagram)
library(eulerr)

set.seed(123)


# since the contrast is different direct validation is not possible, but we can 
#check the overlap of the significant genes in both datasets

# rna - normalised tissue
# prot - direct tumour comparison


# load
# prot
protDE <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/proteomics_results/DE_OCCC_vs_ccRCC_proteomics.csv", col_names = TRUE)
# rna
DE_OCvsRC <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/transcriptomics_results/DE_results_OCCC_vs_ccRCC_ratio.csv")
DE_OCvsGTExOV <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/transcriptomics_results/DE_results_OCCC_vs_GTEx_Ovary_ratio.csv")

# OCCC vs ccRCC
# just sig genes
rna_sig <- DE_OCvsRC$Geneid[DE_OCvsRC$adj.P.Val < 0.05]
prot_sig <- protDE$Geneid[protDE$adj.P.Val < 0.05]


# other
fit <- euler(list(Transcriptomics = rna_sig,
                  Proteomics= prot_sig))

plot(fit,
     fills = list(fill = c("#4393c3", "#d73027"), alpha = 0.5),
     edges = list(col = "black", lwd = 1.5),
     labels = list(fontsize = 11, fontfamily = "sans"),
     quantities = list(fontsize = 10),
     legend = FALSE)


#########################

# MA plot prot
protDE$logFC <- as.numeric(protDE$logFC)
protDE$AveExpr <- as.numeric(protDE$AveExpr)

with(protDE, plot(AveExpr, logFC,
              col = ifelse(adj.P.Val < 0.05, "red", "grey"),
              pch = 16, cex = 0.5))

abline(h = 0, col = "blue")


# MA plot rna
DE_OCvsRC$logFC <- as.numeric(DE_OCvsRC$logFC)
DE_OCvsRC$AveExpr <- as.numeric(DE_OCvsRC$AveExpr)

with(DE_OCvsRC, plot(AveExpr, logFC,
              col = ifelse(adj.P.Val < 0.05, "red", "grey"),
              pch = 16, cex = 0.5))
abline(h = 0, col = "blue")

# save plots with ggplot
ggplot(DE_OCvsRC, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("grey60", "red"),
                     labels = c("Not Significant", "Significant")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  theme_classic() +
  labs(x = "Average log2-expression",
       y = "Log2 FC",
       color = "FDR < 0.05") # adjusted p val

ggplot(protDE, aes(x = AveExpr, y = logFC)) +
  geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c("grey60", "red"),
                     labels = c("Not Significant", "Significant")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  theme_classic() +
  labs(x = "Average log2-expression", y = "Log2 FC",
    color = "FDR < 0.05")

###############################

# multi-omics scatter plot

ggplot(merged, aes(x = logFC.x, y = logFC.y)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    x = "Transcriptomics logFC",
    y = "Proteomics logFC"
  )

cor_val <- cor(merged$logFC.x, merged$logFC.y)
cor_val

pearson_cor  <- cor(merged$logFC.x, merged$logFC.y, method = "pearson")
spearman_cor <- cor(merged$logFC.x, merged$logFC.y, method = "spearman")

summary(merged$logFC.x)
summary(merged$logFC.y)


# save this
cor_val <- cor(merged$logFC.x, merged$logFC.y, 
               method = "spearman", use = "complete.obs")
cor_label <- paste0("ρ = ", round(cor_val, 3)) # spearman rho

ggplot(merged, aes(x = logFC.x, y = logFC.y)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(alpha = 0.4, size = 0.8, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf, label = cor_label, 
           hjust = 1.2, vjust = 1.5,
           size = 3.5) +
  theme_classic() +       
  labs(title = NULL,  x = "Transcriptomics log2 FC",
       y = "Proteomics log2 FC") +
  theme(
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 9),
    plot.margin = margin(10, 15, 10, 10)
  )

