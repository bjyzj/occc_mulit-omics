


library(tidyverse)
library(pheatmap)
library(RColorBrewer)

set.seed(123)

DE_OCvsRC <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/DE_results_OCCC_vs_ccRCC_ratio.csv")
DE_OCvsGTExOV <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/DE_results_OCCC_vs_GTEx_Ovary_ratio.csv")

sample_info <- readRDS("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/sample_info_OCandRC.rds")
expr_mat <- readRDS("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/expression_matrix_occc_renal.rds")

sample_info_ovary_only <- readRDS("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/sample_info_occc_ovary.rds")
expr_mat_ovary_only <- readRDS("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/expression_matrix_occc_ovary.rds")

# OC vs RC
#############
# log2((OCCC / GTEx Ovary) / (ccRCC / GTEx Renal Cortex))
# gtex not included as visuals
# gene selection (50 OCCC-relative + 50 ccRCC-relative)
DE_OCvsRC$adj.P.Val <- pmax(DE_OCvsRC$adj.P.Val, 1e-300)
topN <- 100

sig_up_occc <- DE_OCvsRC %>% 
  filter(adj.P.Val < 0.05 & logFC > 1) %>%
  arrange(adj.P.Val) %>% slice_head(n = topN/2) %>% pull(Geneid)

sig_up_ccrcc <- DE_OCvsRC %>% 
  filter(adj.P.Val < 0.05 & logFC < -1) %>%
  arrange(adj.P.Val) %>% slice_head(n = topN/2) %>% pull(Geneid)

sig_genes <- unique(c(sig_up_occc, sig_up_ccrcc))


tumor_samples <- sample_info$Sample[sample_info$Status == "Tumour"] 
exprDE <- expr_mat[sig_genes, tumor_samples, drop = FALSE]
sample_info_tumor <- sample_info[tumor_samples, ]


occc_cols <- colnames(exprDE)[sample_info_tumor$Tissue == "Ovary"]
ccrcc_cols <- setdiff(colnames(exprDE), occc_cols) 

cat("OCCC samples:", length(occc_cols)) # 845
cat("ccRCC samples:", length(ccrcc_cols)) # 845

expr_scaled <- t(scale(t(exprDE)))
expr_scaled[expr_scaled > 3] <- 3
expr_scaled[expr_scaled < -3] <- -3
expr_scaled[!is.finite(expr_scaled)] <- 0

# row annotation
logFC_vec <- DE_OCvsRC$logFC[match(rownames(expr_scaled), DE_OCvsRC$Geneid)]
row_anno <- data.frame(Direction = ifelse(logFC_vec > 0, "Upregulated in OCCC", "Upregulated in ccRCC"),
                       row.names = rownames(expr_scaled))

annotation_col <- sample_info[match(colnames(expr_scaled), sample_info$Sample),
                              c("Study", "Tissue", "SourceType"), drop = FALSE]

rownames(annotation_col) <- colnames(expr_scaled)
colnames(annotation_col) <- c("Study", "Tissue", "Source Type")


ann_colors <- list(
  Study = c(
    "Bolton et al. (2022)" = "#1f77b4",
    "TCGA-KIRC" = "#B980A7",
    "Expression Atlas (2026)" = "#D54046",
    "GSE160692" = "#9EB258",
    "GSE189553" = "#2E7341",
    "Nagasawa et al. (2019)" = "#D75F22"), 
  Tissue = c(
    "Ovary" = "#17becf",
    "Kidney" = "#CCEDBF"
  ),
  `Source Type` = c(
    "Cell line" = "#E87687",
    "Primary Tumour" = "#F5B35E"
  )
)

row_dist <- as.dist(1 - cor(t(expr_scaled), method = "spearman"))
col_dist <- as.dist(1 - cor(expr_scaled, method = "spearman"))

heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
ann_row_colors <- list(Direction = c("Upregulated in OCCC" = "#B81840","Upregulated in ccRCC" = "#377eb8"))


h1 <- pheatmap(expr_scaled, 
               scale = "none", 
               color = heatmap_colors,
               cluster_rows = TRUE, # hierarchical clustering for rows
               cluster_cols = TRUE, # hierarchical clustering for columns
               clustering_distance_rows = row_dist,
               clustering_distance_cols = col_dist,
               clustering_method = "complete", 
               annotation_col = annotation_col,
               annotation_row = row_anno,
               annotation_colors = c(ann_colors, ann_row_colors),
               show_rownames = FALSE,
               show_colnames = FALSE,
               border_color = NA)
print(h1)
# heatmap_OCvsRC saved
#############################

# top 50 OCCC-relative genes (logFC > 1) with adj.P.Val < 0.05
DE_OCvsGTExOV$adj.P.Val <- pmax(DE_OCvsGTExOV$adj.P.Val, 1e-300)
topN <- 50

sig_up_occc <- DE_OCvsGTExOV %>% filter(adj.P.Val < 0.05 & logFC > 1) %>%
  arrange(adj.P.Val) %>% slice_head(n = topN/2) %>% pull(Geneid)

sig_up_ovary <- DE_OCvsGTExOV %>% filter(adj.P.Val < 0.05 & logFC < -1) %>%
  arrange(adj.P.Val) %>% slice_head(n = topN/2) %>% pull(Geneid)

sig_genes <- unique(c(sig_up_occc, sig_up_ovary))

exprDE_occc <- expr_mat_ovary_only[sig_genes, , drop = FALSE]

# scale rows 
exprDE_scaled <- t(scale(t(exprDE_occc)))
exprDE_scaled[exprDE_scaled > 3] <- 3
exprDE_scaled[exprDE_scaled < -3] <- -3
exprDE_scaled[!is.finite(exprDE_scaled)] <- 0


# row annotation
logFC_vec <- DE_OCvsGTExOV$logFC[match(rownames(exprDE_scaled), DE_OCvsGTExOV$Geneid)]
annotation_row <- data.frame(Direction = ifelse(logFC_vec > 0, "Upregulated in OCCC", "Upregulated in Normal Ovary"),
  row.names = rownames(exprDE_scaled))

# column annotation
annotation_col <- sample_info_ovary_only[match(colnames(exprDE_scaled), sample_info_ovary_only$Sample),
                                         c("Study", "SourceType"), drop = FALSE]
rownames(annotation_col) <- colnames(exprDE_scaled)
colnames(annotation_col) <- c("Study", "Source Type")

# rename GTEx_Ovary to Normal Ovary
annotation_col$Study <- as.character(annotation_col$Study)  # chr
annotation_col$Study[annotation_col$Study == "GTEx_Ovary"] <- "Normal Ovary"
annotation_col$Study <- factor(annotation_col$Study)         # factor

ann_colors <- list(
  Study = c(
    "Bolton et al. (2022)" = "#1f77b4",
    "Normal Ovary" = "#C4D8F3",
    "Expression Atlas (2026)" = "#D54046",
    "GSE160692" = "#9EB258",
    "GSE189553" = "#2E7341",
    "Nagasawa et al. (2019)" = "#D75F22"
  ),
  `Source Type` = c(
    "Cell line" = "#E87687",
    "Primary Tumour" = "#F5B35E"
  )
)


ann_row_colors <- list(Direction = c("Upregulated in OCCC" = "#B81840","Upregulated in Normal Ovary" = "#377eb8"))

dist_rows <- as.dist(1 - cor(t(exprDE_scaled), method = "spearman"))
dist_cols <- as.dist(1 - cor(exprDE_scaled, method = "spearman"))

heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

h2 <- pheatmap(exprDE_scaled, 
               scale = "none",
               color = heatmap_colors, 
               cluster_rows = TRUE,
               cluster_cols = TRUE,  
               clustering_distance_rows = dist_rows,
               clustering_distance_cols = dist_cols,
               clustering_method = "complete",  
               annotation_col = annotation_col,
               annotation_row = annotation_row,
               annotation_colors = c(ann_colors, ann_row_colors),
               show_rownames = FALSE,  
               show_colnames = FALSE,
               border_color = NA)
h2

# no need to dd Group as Study already shows the gtex and occc differnces

# default: clustering_distance_rows = "euclidean": so Euclidean distance + Ward.D2 hierarchical clustering
# i dont need clustering_distance_cols = "correlation" if I set cluster_cols = FALSE, it will cluster by expression similarity not by metadata. 
#If I set cluster_cols = TRUE, it will cluster by metadata and not by expression similarity. So I need to set cluster_cols = FALSE to cluster by expression similarity.


##############################














