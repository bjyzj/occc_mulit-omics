

#install.packages("devtools")
#devtools::install_github("dviraran/xCell")

library(tidyverse)
library(ggplot2)
library(xCell)
library(ggpubr)
library(RColorBrewer)
library(writexl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(GSEABase)
library(pheatmap)

set.seed(123)
# load
#OCCC
CCOC_df_mRNA <- read_tsv("/Users/beyzaerkal/Desktop/occc_multi-omics/inputs/CCOC_RNAseq_mRNA_TPM.tsv")

ea_tpm_mRNA <- read_tsv("/Users/beyzaerkal/Desktop/occc_multi-omics/inputs/E-MTAB_RNAseq_mRNA_TPM.tsv")

tpm_GSE160692_mRNA <- read_tsv("/Users/beyzaerkal/Desktop/occc_multi-omics/inputs/GSE160692_RNA_seq_mRNA_TPM.tsv")

tpm_GSE189553_mRNA <- read_tsv("/Users/beyzaerkal/Desktop/occc_multi-omics/inputs/GSE189553_RNAseq_mRNA_TPM.tsv")

tpm_SD1_OC_mRNA <- read_tsv("/Users/beyzaerkal/Desktop/occc_multi-omics/inputs/SD1_RNAseq_mRNA_TPM.tsv")

#ccrCC
tpm_ccRCC_mRNA <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/input_ccRCC/ccRCC_RNAseq_mRNA_TPM.csv")

# hallmark 
hallmark <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/h.all.v2026.1.Hs.symbols.gmt")

################################


# put genes in rownames
CCOC_df_mRNA <- as.data.frame(CCOC_df_mRNA)
rownames(CCOC_df_mRNA) <- CCOC_df_mRNA$Geneid
CCOC_df_mRNA$Geneid <- NULL

ea_tpm_mRNA <- as.data.frame(ea_tpm_mRNA)
rownames(ea_tpm_mRNA) <- ea_tpm_mRNA$Geneid
ea_tpm_mRNA$Geneid <- NULL

tpm_GSE160692_mRNA <- as.data.frame(tpm_GSE160692_mRNA)
rownames(tpm_GSE160692_mRNA) <- tpm_GSE160692_mRNA$Geneid
tpm_GSE160692_mRNA$Geneid <- NULL

tpm_GSE189553_mRNA <- as.data.frame(tpm_GSE189553_mRNA)
rownames(tpm_GSE189553_mRNA) <- tpm_GSE189553_mRNA$Geneid
tpm_GSE189553_mRNA$Geneid <- NULL

tpm_SD1_OC_mRNA <- as.data.frame(tpm_SD1_OC_mRNA)
rownames(tpm_SD1_OC_mRNA) <- tpm_SD1_OC_mRNA$Geneid
tpm_SD1_OC_mRNA$Geneid <- NULL

tpm_ccRCC_mRNA <- as.data.frame(tpm_ccRCC_mRNA)
rownames(tpm_ccRCC_mRNA) <- tpm_ccRCC_mRNA$Geneid
tpm_ccRCC_mRNA$Geneid <- NULL


# all together
common_occc_genes <- Reduce(intersect, list(
  rownames(CCOC_df_mRNA),
  rownames(ea_tpm_mRNA),
  rownames(tpm_GSE160692_mRNA),
  rownames(tpm_GSE189553_mRNA),
  rownames(tpm_SD1_OC_mRNA)
))

CCOC_df_mRNA <- CCOC_df_mRNA[common_occc_genes, , drop = FALSE]
ea_tpm_mRNA <- ea_tpm_mRNA[common_occc_genes, , drop = FALSE]
tpm_GSE160692_mRNA <- tpm_GSE160692_mRNA[common_occc_genes, , drop = FALSE]
tpm_GSE189553_mRNA <- tpm_GSE189553_mRNA[common_occc_genes, , drop = FALSE]
tpm_SD1_OC_mRNA <- tpm_SD1_OC_mRNA[common_occc_genes, , drop = FALSE]

occc_all <- cbind(
  CCOC_df_mRNA,
  ea_tpm_mRNA,
  tpm_GSE160692_mRNA,
  tpm_GSE189553_mRNA,
  tpm_SD1_OC_mRNA
)

# xCell

exprMatrix <- as.matrix(occc_all)

head(rownames(exprMatrix))

results <- xCellAnalysis(exprMatrix)

dim(results) # 67 235



# sample info for OCCC only
samples <- colnames(occc_all)

ea_samples <- c("ES-2", "JHOC-5", "OVISE", "OVMANA", "OVTOKO", "RMG-I", "TOV-21G")
study <- ifelse(samples %in% ea_samples, "Expression Atlas (2026)",
                ifelse(samples %in% c("C1", "C2", "C3", "C4", "C5", "C6"), "Nagasawa et al. (2019)",
                       ifelse(grepl("^CCC_", samples), "GSE189553",
                              ifelse(grepl("^OVA[0-9]+", samples), "GSE160692",
                                     ifelse(grepl("IGO_07456", samples), "Bolton et al. (2022)",
                                            "Unknown")))))

source_type <- ifelse(study == "Expression Atlas (2026)", "Cell line", "Primary Tumour")

# metadata combined
sample_info <- data.frame(
  Sample = samples,
  Study = factor(study),
  SourceType = factor(source_type)
)


table(sample_info$Study)
sum(is.na(sample_info$Study))

# long format for the heatmap
df_long <- results %>% as.data.frame() %>% rownames_to_column("CellType") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "Score")

unique(df_long$CellType)

ggplot(df_long, aes(x = Sample, y = CellType, fill = Score)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7)) + labs(x = "Samples", y = "Cell Types")


# boxplot setup
df_long <- as.data.frame(results) %>%
  rownames_to_column("CellType") %>%
  pivot_longer(-CellType,
               names_to = "Sample",
               values_to = "Score") %>%
  left_join(sample_info, by = "Sample")

# boxplot
ggplot(df_long %>% filter(CellType == "CD8+ T-cells"),
       aes(x = Study, y = Score, fill = Study)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "CD8 T cell infiltration (xCell)",
       x = "Group",
       y = "Enrichment Score")

ggplot(df_long %>% filter(CellType == "Macrophages"),
       aes(x = Study, y = Score, fill = Study)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "Macrophages",
       x = "Group",
       y = "Enrichment Score")


ggplot(df_long %>% filter(CellType == "Macrophages M2"),
       aes(x = Study, y = Score, fill = Study)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "Macrophages M2",
       x = "Group",
       y = "Enrichment Score")

ggplot(df_long %>% filter(CellType == "StromaScore"),
       aes(x = Study, y = Score, fill = Study)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "StromaScore",
       x = "Group",
       y = "Enrichment Score")


ggplot(df_long %>% filter(CellType == "Fibroblasts"),
       aes(x = Study, y = Score, fill = Study)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "Fibroblasts",
       x = "Group",
       y = "Enrichment Score")

ggplot(df_long %>% filter(CellType == "Fibroblasts"),
       aes(x = Study, y = Score, fill = Study)) +
  geom_boxplot(outlier.size = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(title = "Fibroblasts",
       x = "Group",
       y = "Enrichment Score")


# Is there significant variability across studies
kruskal.test(Score ~ Study, data = df_long %>% filter(CellType == "CD8+ T-cells"))
kruskal.test(Score ~ Study, data = df_long %>% filter(CellType == "CD4+ T-cells"))
kruskal.test(Score ~ Study, data = df_long %>% filter(CellType == "Macrophages"))
kruskal.test(Score ~ Study, data = df_long %>% filter(CellType == "Macrophages M1"))
kruskal.test(Score ~ Study, data = df_long %>% filter(CellType == "Macrophages M2"))
kruskal.test(Score ~ Study, data = df_long %>% filter(CellType == "StromaScore"))
kruskal.test(Score ~ Study, data = df_long %>% filter(CellType == "cDC")) # dendritic cells
kruskal.test(Score ~ Study, data = df_long %>% filter(CellType == "NK cells"))
kruskal.test(Score ~ Study, data = df_long %>% filter(CellType == "Tregs"))

pairwise.wilcox.test(df_long$Score, df_long$Study, p.adjust.method = "BH")







######################
# ccRCC
######################

exprMatrix2 <- as.matrix(tpm_ccRCC_mRNA)

head(rownames(exprMatrix2))

ccrcc_results <- xCellAnalysis(exprMatrix2)

dim(ccrcc_results) # 67 610

# heatmap
df_long_ccrcc <- ccrcc_results %>% as.data.frame() %>% rownames_to_column("CellType") %>%
  pivot_longer(-CellType, names_to = "Sample", values_to = "Score")

ggplot(df_long_ccrcc, aes(x = Sample, y = CellType, fill = Score)) + geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7)) + labs(x = "Samples", y = "Cell Types")



# boxplots



################
# combine
################

df_long_ccrcc$TumourType <- "ccRCC"
df_long$TumourType <- "OCCC"

df_long_ccrcc$Study <- "ccRCC_dataset"  

common_cols <- intersect(names(df_long), names(df_long_ccrcc))
df_long <- df_long[, common_cols]
df_long_ccrcc <- df_long_ccrcc[, common_cols]
df_combined <- bind_rows(df_long, df_long_ccrcc)

cell_types <- c(
  "CD8+ T-cells",
  "CD4+ T-cells",
  "Tregs",
  "NK cells",
  "Macrophages M1",
  "Macrophages M2",
  "StromaScore",
  "ImmuneScore",
  "MicroenvironmentScore",
  "cDC"
)

df_plot <- df_combined %>%
  filter(CellType %in% cell_types)

ggplot(df_plot, aes(x = TumourType, y = Score, fill = TumourType)) +
  geom_boxplot() +
  facet_wrap(~ CellType, scales = "free_y") +
  theme_bw()



# test
p_values <- sapply(cell_types, function(ct) {
  wilcox.test(Score ~ TumourType,
              data = df_combined %>% filter(CellType == ct)
  )$p.value
})

p_adj <- p.adjust(p_values, method = "BH")

results <- data.frame(
  CellType = cell_types,
  p_value = p_values,
  p_adj = p_adj
)

results

#write_xlsx(results, "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/immune_cell_infiltration_results.xlsx")

##############
# occc only biology
##############
df_occ <- df_combined %>% filter(TumourType == "OCCC")

df_occ_mat <- df_occ %>% select(Sample, CellType, Score) %>%
  pivot_wider(names_from = CellType, values_from = Score)



# PCA OCCC + ccRCC
mat <- df_combined %>% dplyr::select(Sample, CellType, Score) %>%
  pivot_wider(names_from = CellType, values_from = Score)

pca <- prcomp(scale(mat[,-1]), center = TRUE)

pca_df <- data.frame(
  Sample = mat$Sample,
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

pca_df <- left_join(pca_df,
                    df_combined %>% distinct(Sample, TumourType),
                    by = "Sample")


pve <- summary(prcomp(scale(mat[,-1])))$importance[2, 1:2] * 100

ggplot(pca_df, aes(PC1, PC2, color = TumourType)) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    color = NA
  ) +
  geom_point(size = 1, alpha = 0.85) +
  
  scale_color_manual(values = c(
    "OCCC" = "#4DBBD5",
    "ccRCC" = "#E64B35"
  )) +
  
  labs(
    x = paste0("PC1 (", round(pve[1], 1), "% variance)"),
    y = paste0("PC2 (", round(pve[2], 1), "% variance)"),
    color = "Tumour Type"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )


##################
# pathway level connection
##################

# OCCC 
exprMatrix_log <- log2(exprMatrix + 1)

# make hallmake list ready format
hallmark_list <- split(hallmark$gene, hallmark$term)
length(hallmark_list)
head(names(hallmark_list))
head(hallmark_list[[1]])


param <- GSVA::ssgseaParam(
  exprData = exprMatrix_log,
  geneSets = hallmark_list,
  minSize = 10,
  maxSize = 500,
  alpha = 0.25,
  normalize = TRUE
)

ssgsea_scores <- gsva(param)


df_ssgsea <- as.data.frame(ssgsea_scores) %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(-Pathway, names_to = "Sample", values_to = "Score")

# heatmap 
mat <- ssgsea_scores
rownames(mat) <- rownames(ssgsea_scores)

pheatmap(mat,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_colnames = FALSE)


pathway_var <- apply(ssgsea_scores, 1, var)
top_paths <- names(sort(pathway_var, decreasing = TRUE))[1:10]
pheatmap(ssgsea_scores[top_paths, ],
         scale = "row")


# correlation between immune cell infiltration and pathway scores
immune_scores <- df_long %>% filter(CellType %in% c("CD8+ T-cells", "Macrophages M1", "Macrophages M2", "Tregs", 
                                                    "StromaScore")) %>% dplyr::select(Sample, CellType, Score) %>% 
  pivot_wider(names_from = CellType, values_from = Score)

pathway_scores <- df_ssgsea %>%
  dplyr::select(Sample, Pathway, Score) %>%
  pivot_wider(names_from = Pathway, values_from = Score)


cor_imm_path <- cor(immune_scores[, -1],
                    pathway_scores[, -1],
                    method = "spearman",
                    use = "pairwise.complete.obs")

pheatmap(cor_imm_path,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")

# size 9 to 5.50
