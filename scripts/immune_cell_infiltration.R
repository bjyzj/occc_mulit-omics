

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
library(circlize)
library(reshape2)
library(ComplexHeatmap)

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
colnames(ssgsea_scores)


ssgsea_scores_df <- as.data.frame(ssgsea_scores)
ssgsea_scores_df$Pathways <- rownames(ssgsea_scores_df)
ssgsea_scores_df <- ssgsea_scores_df[, c("Pathways", setdiff(colnames(ssgsea_scores_df), "Pathways"))]
rownames(ssgsea_scores_df) <- NULL
#write.csv(ssgsea_scores_df, "/Users/beyzaerkal/Desktop/occc_multi-omics/processed/ssgsea_scores.csv", row.names = FALSE)

df_ssgsea <- as.data.frame(ssgsea_scores) %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(-Pathway, names_to = "Sample", values_to = "Score")

# heatmap 
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



##########################
# circo plot

df_long <- results %>%
  as.data.frame() %>%
  rownames_to_column("CellType") %>%
  pivot_longer(-CellType,
               names_to = "Sample",
               values_to = "Score") %>%
  left_join(sample_info, by = "Sample")


cell_study_mat <- df_long %>%
  group_by(CellType, Study) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = "drop")


circos_links <- cell_study_mat %>%
  transmute(
    from = CellType,
    to = Study,
    value = mean_score
  )

threshold <- quantile(circos_links$value, 0.90, na.rm = TRUE)

circos_links <- circos_links %>%
  filter(value >= threshold)

pathway_scores <- df_ssgsea

cell_scores <- df_long %>%
  dplyr::select(CellType, Sample, Score)

pathway_scores <- df_ssgsea %>%
  dplyr::rename(PathwayScore = Score)


cell_mat <- df_long %>%
  group_by(CellType, Sample) %>%
  summarise(cell_score = mean(Score), .groups = "drop") %>%
  pivot_wider(names_from = CellType, values_from = cell_score)

####
hallmark_list <- split(hallmark$gene, hallmark$term)

gene_path_df <- bind_rows(lapply(names(hallmark_list), function(pw) {
  data.frame(
    Gene = hallmark_list[[pw]],
    Pathway = pw,
    value = 1
  )
}))

genes_keep <- rownames(exprMatrix)

gene_path_df <- gene_path_df %>%
  filter(Gene %in% genes_keep) 

gene_path_df_small <- gene_path_df %>%
  group_by(Pathway) %>%
  slice_head(n = 1) %>%  
  ungroup()

circos.clear()

chordDiagram(
  x = gene_path_df_small,
  transparency = 0.5,
  annotationTrack = "grid"
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1] + 0.1,
      CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.5
    )
  },
  bg.border = NA
)






# pathways simialrities only with jaccard index - not sure if I need it
hallmark_list <- split(hallmark$gene, hallmark$term)

pathways <- names(hallmark_list)

overlap_mat <- outer(pathways, pathways, Vectorize(function(p1, p2) {
  g1 <- hallmark_list[[p1]]
  g2 <- hallmark_list[[p2]]
  
  length(intersect(g1, g2)) / length(union(g1, g2))
}))

rownames(overlap_mat) <- pathways
colnames(overlap_mat) <- pathways


overlap_df <- melt(overlap_mat)
colnames(overlap_df) <- c("from", "to", "value")

overlap_df <- overlap_df %>%
  filter(from != to)

threshold <- quantile(overlap_df$value, 0.95, na.rm = TRUE)

overlap_df_filtered <- overlap_df %>%
  filter(value >= threshold)

circos.clear()

chordDiagram(
  overlap_df_filtered,
  transparency = 0.4,
  annotationTrack = "grid"
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1] + 0.1,
      CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      cex = 0.5
    )
  },
  bg.border = NA
)








####################
##########################
# Circos: Genes --> Pathways --> Immune Cells
##########################

# ── Step 1: build aligned matrices (OCCC samples only) ───────────────────────
# df_long at this point is OCCC-only (it was rebuilt from results + sample_info)
# ssgsea was also run on exprMatrix (OCCC only) so samples should match

cell_mat <- df_long %>%
  filter(CellType %in% c("CD8+ T-cells", "CD4+ T-cells", "Tregs",
                         "NK cells", "Macrophages M1", "Macrophages M2",
                         "StromaScore")) %>%
  group_by(Sample, CellType) %>%
  summarise(Score = mean(Score), .groups = "drop") %>%
  pivot_wider(names_from = CellType, values_from = Score)

immune_cells_of_interest <- c(
  "CD8+ T-cells", "CD4+ T-cells", "Tregs",
  "NK cells", "Macrophages M1", "Macrophages M2",
  "StromaScore"
)

all_immune <- immune_cells_of_interest

immune_cols <- setNames(
  brewer.pal(length(all_immune), "Set2"),
  all_immune
)


path_mat <- df_ssgsea %>%
  pivot_wider(names_from = Pathway, values_from = Score)

# align rows by Sample — CRITICAL, otherwise cor() is garbage
shared_samples <- intersect(cell_mat$Sample, path_mat$Sample)
stopifnot(length(shared_samples) > 5)   # sanity check

cell_mat  <- cell_mat[match(shared_samples, cell_mat$Sample), ]
path_mat  <- path_mat[match(shared_samples, path_mat$Sample), ]

# ── Step 2: immune cell × pathway correlation ─────────────────────────────────
# result: rows = immune cell types, cols = hallmark pathways
immune_path_cor <- cor(
  cell_mat[ , -1],   # drop Sample column
  path_mat[ , -1],
  method = "spearman",
  use    = "pairwise.complete.obs"
)

# ── Step 3: select top pathways per immune cell (|r| >= 0.4) ─────────────────
cor_long <- immune_path_cor %>%
  as.data.frame() %>%
  rownames_to_column("CellType") %>%
  pivot_longer(-CellType, names_to = "Pathway", values_to = "r") %>%
  filter(abs(r) >= 0.4)                         # tune threshold as needed

# ── Step 4: map pathways to their top genes via ssGSEA contribution ───────────
# use hallmark_list (already built earlier) 
# for each retained pathway, take the top 5 genes by mean expression in OCCC
mean_expr <- rowMeans(exprMatrix_log)            # named vector, gene → mean log2TPM

gene_path_links <- bind_rows(lapply(unique(cor_long$Pathway), function(pw) {
  genes_in_pw <- intersect(hallmark_list[[pw]], names(mean_expr))
  if (length(genes_in_pw) == 0) return(NULL)
  top_genes <- sort(mean_expr[genes_in_pw], decreasing = TRUE)
  top_genes <- head(top_genes, 5)               # top 5 expressed genes per pathway
  data.frame(
    from  = names(top_genes),                   # Gene
    to    = pw,                                  # Pathway
    value = as.numeric(top_genes),              # chord width = mean expression
    stringsAsFactors = FALSE
  )
}))

# pathway → immune cell links (use |r| as chord width)
path_immune_links <- cor_long %>%
  transmute(
    from  = Pathway,
    to    = CellType,
    value = abs(r)
  )

# combined edge table for chordDiagram
chord_df <- bind_rows(gene_path_links, path_immune_links)

# ── Step 5: sector colours ────────────────────────────────────────────────────
all_genes      <- unique(gene_path_links$from)
all_pathways   <- unique(c(gene_path_links$to, path_immune_links$from))
all_immune     <- unique(path_immune_links$to)

n_pw <- length(all_pathways)
n_im <- length(all_immune)
n_ge <- length(all_genes)

# colour palettes per sector group
gene_cols <- setNames(
  colorRampPalette(brewer.pal(9, "Set1"))(n_ge),
  all_genes
)

# pathways: warm palette
pathway_cols <- setNames(
  colorRampPalette(brewer.pal(9, "YlOrRd"))(n_pw),
  all_pathways
)

immune_cols  <- setNames(
  brewer.pal(max(3, n_im), "Set2")[seq_len(n_im)], all_immune)   # qualitative for immune

grid_col <- c(gene_cols, pathway_cols, immune_cols)

# ── Step 6: which sectors to label (skip genes to avoid clutter) ─────────────
label_sectors <- c(all_genes, all_pathways, all_immune)
# ── Step 7: signed colour lookup for pathway sector annotation ───────────────
# average correlation of each pathway across all immune cells
pathway_avg_cor <- cor_long %>%
  group_by(Pathway) %>%
  summarise(avg_r = mean(r), .groups = "drop")

col_fun <- colorRamp2(c(-0.6, 0, 0.6), c("#4575B4", "white", "#D73027"))

# ── Step 8: draw ──────────────────────────────────────────────────────────────
circos.clear()
circos.par(
  gap.after = c(
    rep(1, n_ge  - 1), 8,   # gap after gene block
    rep(1, n_pw  - 1), 8,   # gap after pathway block
    rep(1, n_im  - 1), 8    # gap after immune block
  ),
  start.degree = 90
)

chordDiagram(
  chord_df,
  grid.col          = grid_col,
  transparency      = 0.45,
  annotationTrack   = "grid",
  preAllocateTracks = list(list(track.height = mm_h(5)))  # ONE pre-allocated track
)

# ── Step 9: single panel.fun — rects + labels in ONE call (fixes your bug) ───
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    
    sector <- CELL_META$sector.index
    
    # colour rect for pathway sectors
    if (sector %in% all_pathways) {
      avg_r    <- pathway_avg_cor$avg_r[pathway_avg_cor$Pathway == sector]
      rect_col <- if (length(avg_r) > 0) col_fun(avg_r) else "grey80"
      circos.rect(
        CELL_META$xlim[1], CELL_META$ylim[1],
        CELL_META$xlim[2], CELL_META$ylim[2],
        col    = rect_col,
        border = NA
      )
    }
    
    # labels for ALL sectors (genes + pathways + immune)
    if (sector %in% label_sectors) {
      
      # clean up pathway names; leave gene/immune names as-is
      label <- if (sector %in% all_pathways) {
        gsub("_", " ", gsub("HALLMARK_", "", sector))
      } else {
        sector
      }
      
      circos.text(
        CELL_META$xcenter,
        CELL_META$ylim[1] + 0.1,
        label,
        facing     = "clockwise",
        niceFacing = TRUE,
        adj        = c(0, 0.5),
        cex        = if (sector %in% all_genes) 0.4 else 0.45  # genes slightly smaller
      )
    }
  },
  bg.border = NA
)
# ── Step 10: legend ───────────────────────────────────────────────────────────
lgd_immune <- Legend(
  labels    = all_immune,
  type      = "point",
  pch       = 16,
  legend_gp = gpar(col = immune_cols),
  title     = "Immune Cell Type",
  title_gp  = gpar(fontsize = 9, fontface = "bold"),
  labels_gp = gpar(fontsize = 8)
)

lgd_cor <- Legend(
  col_fun   = col_fun,
  title     = "Avg. Spearman r\n(pathway–immune)",
  title_gp  = gpar(fontsize = 9, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  direction = "horizontal",
  at        = c(-0.6, -0.3, 0, 0.3, 0.6)
)

draw(packLegend(lgd_immune, lgd_cor, direction = "vertical"),
     x = unit(1, "npc") - unit(2, "mm"),
     y = unit(0.5, "npc"),
     just = c("right", "center"))

circos.clear()
