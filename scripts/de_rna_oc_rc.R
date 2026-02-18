
#################
# RNA-seq only
#################

library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(pheatmap)
library(circlize)
library(limma)
library(readxl)
library(ggrepel)
library(EnhancedVolcano)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(RColorBrewer)
library(igraph)
library(ggraph)
library(cowplot)
library(grid)
library(gridExtra)
library(patchwork)
library(ggplotify)
library(pdftools)
library(magick)
library(qpdf)
library(ggtext)

set.seed(222)

# OVARY
CCOC_df_mRNA <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/CCOC_RNAseq_mRNA_TPM.csv") # 200

ea_tpm_mRNA <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/E-MTAB_RNAseq_mRNA_TPM.csv") # 7

tpm_GSE160692_mRNA <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE160692_RNA_seq_mRNA_TPM.csv") # 11

tpm_GSE189553_mRNA <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE189553_RNAseq_mRNA_TPM.csv") # 11

tpm_113_mRNA <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/tpm113_RNAseq_mRNA_TPM.csv") # 105

tpm_SD1_OC_mRNA <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/SD1_RNAseq_mRNA_TPM.csv") # 6

# RENAL
tpm_ccRCC <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/input_ccRCC/ccRCC_RNAseq_mRNA_TPM.csv")

# GTEx
ovary_gtex <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/raw_gtex02/gene_tpm_v11_ovary.gct", skip = 2)
renal_gtex <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/raw_gtex02/gene_tpm_v11_kidney_cortex.gct", skip = 2)


######################
tpm_SD1_OC_mRNA <- as.data.frame(tpm_SD1_OC_mRNA)
rownames(tpm_SD1_OC_mRNA) <- make.unique(as.character(tpm_SD1_OC_mRNA$Geneid)) # set row names from colnames
#tpm_SD1_OC_mRNA$Geneid <- NULL # delete col

tpm_113_mRNA <- as.data.frame(tpm_113_mRNA)
rownames(tpm_113_mRNA) <- make.unique(as.character(tpm_113_mRNA$Geneid))
#tpm_113_mRNA$Geneid <- NULL

tpm_GSE189553_mRNA <- as.data.frame(tpm_GSE189553_mRNA)
rownames(tpm_GSE189553_mRNA) <- make.unique(as.character(tpm_GSE189553_mRNA$Geneid))
#tpm_GSE189553_mRNA$Geneid <- NULL

ea_tpm_mRNA <- as.data.frame(ea_tpm_mRNA)
rownames(ea_tpm_mRNA) <- make.unique(as.character(ea_tpm_mRNA$Geneid))
#ea_tpm_mRNA$Geneid <- NULL

tpm_GSE160692_mRNA <- as.data.frame(tpm_GSE160692_mRNA)
rownames(tpm_GSE160692_mRNA) <- make.unique(as.character(tpm_GSE160692_mRNA$Geneid))
#tpm_GSE160692_mRNA$Geneid <- NULL

CCOC_df_mRNA <- as.data.frame(CCOC_df_mRNA)
rownames(CCOC_df_mRNA) <- make.unique(as.character(CCOC_df_mRNA$Geneid))
#CCOC_df_mRNA$Geneid <- NULL
##########################
# GTEx clean - ovary

ovary_gtex2 <- ovary_gtex %>%
  dplyr::select(-1) %>%
  dplyr::rename(Geneid = Description) %>%
  dplyr::mutate(Geneid = sub("\\..*", "", Geneid))

#listFilters(ensembl)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# HGNC symbols
bm_results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                    filters = "ensembl_gene_id",
                    values = ovary_gtex2$Geneid,   
                    mart = ensembl)

ovary_gtex_annot <- ovary_gtex2 %>%
  dplyr::left_join(bm_results, by = c("Geneid" = "ensembl_gene_id")) %>%
  dplyr::mutate(Geneid = ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "", hgnc_symbol, Geneid)) %>%
  dplyr::select(-hgnc_symbol)


protein_coding_genes <- getBM(attributes = c("hgnc_symbol"),
                              filters = "biotype",
                              values = "protein_coding",
                              mart = ensembl)
# ENSG IDs without a mapped symbol are removed - only protein-coding genes that have an HGNC symbol
ovary_gtex_filtered <- ovary_gtex_annot %>% filter(Geneid %in% protein_coding_genes$hgnc_symbol)

sum(grepl("^ENSG", ovary_gtex2$Geneid))
sum(!grepl("^ENSG", ovary_gtex2$Geneid))

sum(grepl("^ENSG", ovary_gtex_filtered$Geneid))
sum(!grepl("^ENSG", ovary_gtex_filtered$Geneid)) # 19072 genes

dim(ovary_gtex_filtered) # 19276  195

#####################################
# GTEx clean - renal

renal_gtex2 <- renal_gtex %>%
  dplyr::select(-1) %>%
  dplyr::rename(Geneid = Description) %>%
  dplyr::mutate(Geneid = sub("\\..*", "", Geneid))

bm_results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                    filters = "ensembl_gene_id",
                    values = renal_gtex2$Geneid,   
                    mart = ensembl)

renal_gtex_annot <- renal_gtex2 %>%
  dplyr::left_join(bm_results, by = c("Geneid" = "ensembl_gene_id")) %>%
  dplyr::mutate(Geneid = ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "", hgnc_symbol, Geneid)) %>%
  dplyr::select(-hgnc_symbol)


protein_coding_genes <- getBM(attributes = c("hgnc_symbol"),
                              filters = "biotype",
                              values = "protein_coding",
                              mart = ensembl)

renal_gtex_filtered <- renal_gtex_annot %>% filter(Geneid %in% protein_coding_genes$hgnc_symbol)

sum(grepl("^ENSG", renal_gtex2$Geneid))
sum(!grepl("^ENSG", renal_gtex2$Geneid))

sum(grepl("^ENSG", renal_gtex_filtered$Geneid))
sum(!grepl("^ENSG", renal_gtex_filtered$Geneid))

dim(renal_gtex_filtered) # 19276   106

###################################
# combine OCCC studies
clean_genes <- function(df){
  df <- as.data.frame(df)
  df$Geneid <- toupper(trimws(df$Geneid)) # remove spaces and uppercase all
  df <- df[!duplicated(df$Geneid), ] # remove duplicates, keep first
  rownames(df) <- df$Geneid
  df$Geneid <- NULL
  return(df)
}

CCOC_df_mRNA <- clean_genes(CCOC_df_mRNA)
ea_tpm_mRNA <- clean_genes(ea_tpm_mRNA)
tpm_GSE160692_mRNA <- clean_genes(tpm_GSE160692_mRNA)
tpm_GSE189553_mRNA <- clean_genes(tpm_GSE189553_mRNA)
tpm_113_mRNA <- clean_genes(tpm_113_mRNA) # only 150 Gene symbols decided not to include other 17K
tpm_SD1_OC_mRNA <- clean_genes(tpm_SD1_OC_mRNA)

nrow(CCOC_df_mRNA)
nrow(ea_tpm_mRNA)
nrow(tpm_GSE160692_mRNA)
nrow(tpm_GSE189553_mRNA)
nrow(tpm_SD1_OC_mRNA)
nrow(tpm_113_mRNA)

ncol(CCOC_df_mRNA)
ncol(ea_tpm_mRNA)
ncol(tpm_GSE160692_mRNA)
ncol(tpm_GSE189553_mRNA)
ncol(tpm_SD1_OC_mRNA)
ncol(tpm_113_mRNA)

common_occc_genes <- Reduce(intersect, list(
  rownames(CCOC_df_mRNA),
  rownames(ea_tpm_mRNA),
  rownames(tpm_GSE160692_mRNA),
  rownames(tpm_GSE189553_mRNA),
  rownames(tpm_SD1_OC_mRNA)
))

CCOC_df_mRNA <- CCOC_df_mRNA[common_occc_genes, ]
ea_tpm_mRNA <- ea_tpm_mRNA[common_occc_genes, ]
tpm_GSE160692_mRNA <- tpm_GSE160692_mRNA[common_occc_genes, ]
tpm_GSE189553_mRNA <- tpm_GSE189553_mRNA[common_occc_genes, ]
tpm_SD1_OC_mRNA <- tpm_SD1_OC_mRNA[common_occc_genes, ]

occc_all <- cbind(
  CCOC_df_mRNA,
  ea_tpm_mRNA,
  tpm_GSE160692_mRNA,
  tpm_GSE189553_mRNA,
  tpm_SD1_OC_mRNA
)
######################
anyDuplicated(ovary_gtex_filtered$Geneid)
anyDuplicated(renal_gtex_filtered$Geneid)

# handle duplicate geneid error - take mean of the duplicates
ovary_gtex_filtered <- ovary_gtex_filtered %>%
  group_by(Geneid) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()

anyDuplicated(ovary_gtex_filtered$Geneid)


renal_gtex_filtered <- renal_gtex_filtered %>%
  group_by(Geneid) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup()

anyDuplicated(renal_gtex_filtered$Geneid)
#########################
# geneid col as non numeric clean for gtex 
ovary_gtex_filtered <- as.data.frame(ovary_gtex_filtered)
rownames(ovary_gtex_filtered) <- ovary_gtex_filtered[,1]
ovary_gtex_filtered <- ovary_gtex_filtered[,-1]

renal_gtex_filtered <- as.data.frame(renal_gtex_filtered)
rownames(renal_gtex_filtered) <- renal_gtex_filtered[,1]
renal_gtex_filtered <- renal_gtex_filtered[,-1]

tpm_ccRCC <- as.data.frame(tpm_ccRCC)
rownames(tpm_ccRCC) <- tpm_ccRCC[,1]
tpm_ccRCC <- tpm_ccRCC[,-1]

# log2 transform all datasets (after adding pseudocount of 1)
occc_log <- log2(occc_all + 1)
ovary_gtex_log <- log2(ovary_gtex_filtered + 1)
ccrcc_log <- log2(tpm_ccRCC + 1)
renal_gtex_log <- log2(renal_gtex_filtered + 1)

######################
# expression matrix of both ovary and renal

common_genes <- Reduce(intersect, list(
  rownames(occc_log),
  rownames(ovary_gtex_log),
  rownames(ccrcc_log),
  rownames(renal_gtex_log)
))

occc_log <- occc_log[common_genes, ]
ovary_gtex_log <- ovary_gtex_log[common_genes, ]
ccrcc_log <- ccrcc_log[common_genes, ]
renal_gtex_log <- renal_gtex_log[common_genes, ]
# remove extra col info from gtex
ovary_gtex_log <- ovary_gtex_log[, !colnames(ovary_gtex_log) %in% "gene_biotype"]
renal_gtex_log <- renal_gtex_log[, !colnames(renal_gtex_log) %in% "gene_biotype"]

expr_mat <- cbind(
  occc_log,
  ovary_gtex_log,
  ccrcc_log,
  renal_gtex_log
)

##############
samples <- colnames(expr_mat)

# GTEx ovary and renal sample IDs - for separate labelling

ovary_samples <- colnames(ovary_gtex_log)
renal_samples <- colnames(renal_gtex_log)
ccrcc_samples <- colnames(ccrcc_log)

ea_samples <- c("ES-2", "JHOC-5", "OVISE", "OVMANA", "OVTOKO", "RMG-I", "TOV-21G")

study <- ifelse(samples %in% ea_samples, "Expression Atlas (2026)",
                ifelse(samples %in% c("C1", "C2", "C3", "C4", "C5", "C6"), "Nagasawa et al. (2019)",
                       ifelse(grepl("^CCC_", samples), "GSE189553",
                              ifelse(grepl("^OVA[0-9]+", samples), "GSE160692",
                                     ifelse(grepl("IGO_07456", samples), "Bolton et al. (2022)",
                                            ifelse(samples %in% ccrcc_samples, "ccRCC",
                                                   ifelse(samples %in% ovary_samples, "GTEx_Ovary",
                                                          ifelse(samples %in% renal_samples, "GTEx_Kidney", NA))))))))

tissue <- ifelse(samples %in% ovary_samples | grepl("^ES-|^C[1-6]|^CCC_|^OVA|^IGO", samples), "Ovary", "Kidney")
status <- ifelse(samples %in% ovary_samples | samples %in% renal_samples, "Normal", "Tumor")
source_type <- ifelse(study == "Expression Atlas (2026)", "cell_line", "primary_tumour")

# metadata combined
sample_info <- data.frame(
  Sample = samples,
  Study = factor(study),
  Tissue = factor(tissue),
  Status = factor(status),
  SourceType = factor(source_type)
)
rownames(sample_info) <- sample_info$Sample
samples[duplicated(samples)] # col checks that have duplicates
all(colnames(expr_mat) == rownames(sample_info)) 

# Combined group factor - make categorical data combining metadata into groups
group <- rep(NA, length(samples))
names(group) <- samples
group[samples %in% colnames(occc_log)] <- "OCCC"
group[samples %in% colnames(ovary_gtex_log)] <- "GTEx_Ovary"
group[samples %in% colnames(ccrcc_log)] <- "ccRCC"
group[samples %in% colnames(renal_gtex_log)] <- "GTEx_Kidney"

sample_info$Group <- factor(group, levels = c("GTEx_Ovary", "OCCC", "GTEx_Kidney", "ccRCC"))

table(sample_info$Group)
table(sample_info$Tissue, sample_info$Status)

##############################
# gene checks
# Use top 500 most variable genes 
top_genes <- head(order(apply(expr_mat, 1, sd), decreasing = TRUE), 500)
expr_sub <- expr_mat[top_genes, ]

# Create annotation for columns
annotation_col <- data.frame(
  Group = sample_info$Group,
  Tissue = sample_info$Tissue,
  Status = sample_info$Status
)
rownames(annotation_col) <- sample_info$Sample

pheatmap(expr_sub,
         scale = "row",             
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         clustering_method = "complete")


###############################
# log2((OCCC / GTEx Ovary) / (ccRCC / GTEx Renal Cortex))

design <- model.matrix(~0 + Group, data = sample_info)
colnames(design) <- levels(sample_info$Group)


contrast.matrix <- makeContrasts(
  OCCC_vs_ccRCC_ratio = (OCCC - GTEx_Ovary) - (ccRCC - GTEx_Kidney),
  levels = design
)

fit <- lmFit(expr_mat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

topTable(fit2, coef="OCCC_vs_ccRCC_ratio", number=Inf, adjust.method="BH")


################################
# FILTER LOW EXPRESSED GENES: Keep genes with log2(TPM+1) > 1
N_MIN_SAMPLES <- 6
keep_genes <- rowSums(expr_mat > 1) >= N_MIN_SAMPLES
expr_filtered <- expr_mat[keep_genes, ]

occc_studies <- c("EA","SD1","S113","GSE189553","GSE160692","CCOC")
keep <- sample_info_all$Study %in% c(occc_studies, "ccRCC")

expr_sub <- expr_filtered[, keep]
sample_info_sub <- sample_info_all[keep, ]
# remove any genes that have zero variance after filtering
zero_var <- apply(expr_sub, 1, var) == 0
expr_sub <- expr_sub[!zero_var, ]


bio_group <- factor(ifelse(sample_info_sub$Study %in% occc_studies, "OCCC", "ccRCC")) # check if in the list
batch <- factor(sample_info_sub$Study)
table(bio_group, batch)
design <- model.matrix(~ 0 + bio_group) # in separate col, 


colnames(design)[1:length(levels(bio_group))] <- levels(bio_group)
##################################
#limma
fit <- lmFit(expr_sub, design)
# getting occ and ccr
contrast_matrix <- makeContrasts(OCCC_vs_ccRCC = OCCC - ccRCC, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, coef = "OCCC_vs_ccRCC", number = Inf, adjust.method = "BH")
summary(results$adj.P.Val < 0.05) # 4579   14761

#write.csv(results, file = "/Users/beyzaerkal/Desktop/occc_repo/DE_results_OCCC_vs_ccRCC.csv", row.names = TRUE)
# save sig genes
sig_genes <- results[results$adj.P.Val < 0.05, ] 
#write.csv(sig_genes, file ="/Users/beyzaerkal/Desktop/occc_repo/DE_sig_genes.csv", row.names = TRUE)

up_genes <- results[results$adj.P.Val < 0.05 & results$logFC > 0, ]
down_genes <- results[results$adj.P.Val < 0.05 & results$logFC < 0, ]

#write.table(rownames(up_genes), file ="/Users/beyzaerkal/Desktop/occc_repo/up_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE) 
#write.table(rownames(down_genes), file ="/Users/beyzaerkal/Desktop/occc_repo/down_genes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE) 

################# HEATMAP 1 ""without"" batch correction scaled tgt
topN <- 500
top_genes <- rownames(results[order(results$adj.P.Val), ])[1:topN]

mat <- expr_sub[top_genes, , drop = FALSE]

mat_scaled <- t(scale(t(mat)))
mat_scaled[mat_scaled > 3] <- 3
mat_scaled[mat_scaled < -3] <- -3

study_map <- c(
  CCOC = "OCCRNA1",
  ccRCC = "RCCRNA1",
  EA = "OCCRNA2",
  GSE160692 = "OCCRNA3",
  GSE189553 = "OCCRNA4",
  SD1 = "OCCRNA5",
  S113 = "OCCRNA6")

sample_info_sub$Study_ID <- as.character(sample_info_sub$Study)
sample_info_sub$Study_ID <- ifelse(sample_info_sub$Study_ID %in% names(study_map),
                                study_map[sample_info_sub$Study_ID],
                                sample_info_sub$Study_ID)

sample_info_sub$SourceType <- as.character(sample_info_sub$SourceType)
tumour_type <- as.character(ifelse(sample_info_sub$Study %in% occc_studies, "OCCC", "ccRCC"))


unique_studies <- unique(sample_info_sub$Study_ID)
study_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Set3"))(length(unique_studies)),
  unique_studies
)

bio_colors <- c("OCCC" = "#d73027", "ccRCC" = "#4575b4")
source_colors <- c("cell_line" = "#fdae61", "primary_tumour" = "#313695")

ha_col <- HeatmapAnnotation(
  Study = sample_info_sub$Study_ID,
  TumourType = tumour_type,
  SourceType = sample_info_sub$SourceType,
  col = list(
    Study = study_colors,
    TumourType = bio_colors,
    SourceType = source_colors
  ),
  
  annotation_name_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  annotation_legend_param = list(
    Study = list(title = "Study", ncol = 2,
                 title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                 labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
    TumourType = list(title = "Tumour Type", 
                      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
    SourceType = list(title = "Source Type",
                      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica"))
  )
)


row_anno <- data.frame(Direction = ifelse(results[top_genes, "logFC"] > 0, "Up_in_OCCC", "Up_in_ccRCC"))
rownames(row_anno) <- top_genes
row_colors <- c("Up_in_OCCC" = "#e41a1c", "Up_in_ccRCC" = "#377eb8")
ha_row <- rowAnnotation(Direction = row_anno$Direction, col = list(Direction = row_colors))
# without batch correction - OCCC vs ccRCC Gene Expression
Heatmap(
  mat_scaled,
  name = "z-score",
  top_annotation = ha_col,
  #right_annotation = ha_row,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = function(x) as.dist(1 - cor(t(x), method = "spearman")),
  clustering_distance_columns = function(x) as.dist(1 - cor(x, method = "spearman")),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 7, fontfamily = "Helvetica"),
  row_title = paste0("Top ", topN, " DE genes"),
  row_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  heatmap_legend_param = list(title = "Expression (Z-score)",
                              title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                              labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica"))
)



table(results$logFC > 0) 

table(results[order(results$adj.P.Val), ][1:300, "logFC"] > 0)

table(bio_group, sample_info_sub$Study)


######################################## HEATMAP 2 without batch

# per group sclaed OCCC adn ccRCC. separatluy 
# top 50 OC and 50 RC genes with most sig and z score speartly, set limit for extreme val
topN <- 100
top_up_occc <- results %>% filter(logFC > 0) %>% arrange(adj.P.Val) %>%
  head(topN/2) %>%
  rownames()

top_up_ccrcc <- results %>% filter(logFC < 0) %>% arrange(adj.P.Val) %>%
  head(topN/2) %>%
  rownames()

results[top_up_occc, c("logFC", "adj.P.Val")]
results[top_up_ccrcc, c("logFC", "adj.P.Val")]


top_genes <- c(top_up_occc, top_up_ccrcc)
top_genes <- top_genes[top_genes %in% rownames(expr_sub)] 

mat <- expr_sub[top_genes, , drop = FALSE]

samples_to_plot <- colnames(mat)
rownames(sample_info_sub) <- sample_info_sub$Sample
sample_info_subset <- sample_info_sub[samples_to_plot, , drop = FALSE]
# bio group -> tumourtype
bio_group_subset <- factor(ifelse(sample_info_subset$Study %in% occc_studies, "OCCC", "ccRCC"))
source_type_subset <- factor(sample_info_subset$SourceType)

# row annot - direction

row_anno <- data.frame(Direction = ifelse(results[rownames(mat), "logFC"] > 0, "Up in OCCC", "Up in ccRCC"))
rownames(row_anno) <- rownames(mat)
row_anno$Direction <- factor(row_anno$Direction, levels = c("Up in OCCC", "Up in ccRCC"))
row_colors <- c("Up in OCCC" = "#e41a1c", "Up in ccRCC" = "#377eb8")
ha_row <- rowAnnotation(Direction = row_anno$Direction, col = list(Direction = row_colors))

# scale separtly
mat_scaled <- mat
mat_scaled[, bio_group_subset == "OCCC"] <- t(scale(t(mat[, bio_group_subset == "OCCC"])))
mat_scaled[, bio_group_subset == "ccRCC"] <- t(scale(t(mat[, bio_group_subset == "ccRCC"])))

# cap extremes
mat_scaled[mat_scaled > 3] <- 3
mat_scaled[mat_scaled < -3] <- -3
# remove 
mat_scaled <- mat_scaled[apply(mat_scaled, 1, function(x) all(is.finite(x))), ]

row_anno <- row_anno[rownames(mat_scaled), , drop = FALSE]
ha_row <- rowAnnotation(Direction = row_anno$Direction,
                        col = list(Direction = row_colors), 
                        annotation_name_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                        annotation_legend_param = list(
                          Direction = list(title = "Gene Direction",
                                           title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                                           labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica"))
                        )
)

unique_studies <- unique(sample_info_subset$Study_ID)
study_colors <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(unique_studies)), unique_studies)
bio_colors <- c("OCCC" = "#d73027", "ccRCC" = "#4575b4")
source_colors <- c("cell_line" = "#fdae61", "primary_tumour" = "#313695")

ha_col <- HeatmapAnnotation(
  Study = sample_info_subset$Study_ID,
  TumourType = bio_group_subset,
  SourceType = source_type_subset,
  col = list(
    Study = study_colors,
    TumourType = bio_colors,
    SourceType = source_colors
  ),
  annotation_name_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  annotation_legend_param = list(
    Study = list(title = "Study", ncol = 2, 
                 title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                 labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
    TumourType = list(title = "Tumour Type", 
                      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
    SourceType = list(title = "Source Type",
                      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")
    )
  )
)


col_fun <- colorRamp2(c(-3, 0, 3), c("#377eb8", "white", "#e41a1c"))
# OCCC vs ccRCC expression pattern (without Batch-correction)
h2 <- Heatmap(
  mat_scaled,
  name = "z-score",
  col = col_fun,
  top_annotation = ha_col,
  right_annotation = ha_row,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = function(x) as.dist(1 - cor(t(x), method = "spearman", use = "pairwise.complete.obs")),
  clustering_distance_columns = function(x) as.dist(1 - cor(x, method = "spearman", use = "pairwise.complete.obs")),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 7, fontfamily = "Helvetica"),
  row_title = paste0("Top ", topN, " DE genes"),
  row_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  heatmap_legend_param = list(title = "Expression (Z-score)",
                              title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                              labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica"))
)
print(h2)

###########################################
# with batch correction and repeated
expr_occc <- expr_sub[, sample_info_sub$Study %in% occc_studies]
batch_occc <- factor(sample_info_sub$Study[sample_info_sub$Study %in% occc_studies])

# remove zero-variance genes in OCCC before batch correction
keep_genes_occc <- apply(expr_occc, 1, var) > 0
expr_occc <- expr_occc[keep_genes_occc, ]

expr_occc_corrected <- removeBatchEffect(expr_occc, batch = batch_occc)

# ccRCC expression
expr_ccrcc <- expr_sub[, sample_info_sub$Study == "ccRCC"]

# common genes
common_genes <- intersect(rownames(expr_occc_corrected), rownames(expr_ccrcc))
expr_occc_corrected <- expr_occc_corrected[common_genes, ]
expr_ccrcc <- expr_ccrcc[common_genes, ]

expr_combined <- cbind(expr_occc_corrected, expr_ccrcc)

# remove zero-variance genes fro combined ones
keep_genes_all <- apply(expr_combined, 1, var) > 0
expr_corrected <- expr_combined[keep_genes_all, ]

# batch corrected expression data
#write.csv(expr_corrected, file ="/Users/beyzaerkal/Desktop/occc_repo/Batch_corrected_CC_RC.csv", row.names = TRUE)
####################################
# pca
expr_sub_filtered <- expr_sub[apply(expr_sub, 1, var) > 0, ]

pca_before <- prcomp(t(expr_sub_filtered), scale. = TRUE)
pca_after <- prcomp(t(expr_corrected), scale. = TRUE)

bio_group_char <- ifelse(sample_info_sub$Study %in% occc_studies, "OCCC", "ccRCC")

summary(pca_before)
summary(pca_after)

pca_before_df <- data.frame(pca_before$x[, 1:2], Study = sample_info_sub$Study, TumourType = bio_group_char)
pca_after_df <- data.frame(pca_after$x[, 1:2], Study = sample_info_sub$Study, TumourType = bio_group_char)

ggplot(pca_before_df, aes(PC1, PC2, color = Study, shape = TumourType)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA Before Batch Correction")

ggplot(pca_after_df, aes(PC1, PC2, color = Study, shape = TumourType)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA After Batch Correction")

# set size for  both at the smae time
common_theme <- theme_minimal(base_size = 7) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(size = 7, face = "bold"),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6)
  )


# PCA before batch correction
pca_before_plot <- ggplot(pca_before_df, aes(PC1, PC2, color = TumourType, fill = TumourType)) +
  geom_point(size = 1.2, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, alpha = 0.2, geom = "polygon") + 
  common_theme +
  ggtitle("PCA Before Batch Correction") +
  labs(color = "Tumour Type", fill = "Tumour Type") +
  theme(
    legend.position = "right",
    panel.grid = element_blank())

# PCA after batch correction
pca_after_plot <- ggplot(pca_after_df, aes(PC1, PC2, color = TumourType, fill = TumourType)) +
  geom_point(size = 1.2, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.95, alpha = 0.2, geom = "polygon") + 
  common_theme +
  ggtitle("PCA After Batch Correction") +
  labs(color = "Tumour Type", fill = "Tumour Type") +
  theme(
    legend.position = "right",
    panel.grid = element_blank())
# combine vertically
combined_plot <- plot_grid(pca_before_plot, pca_after_plot, ncol = 1, align = "v", labels = "A")
#ggsave("/Users/beyzaerkal/Desktop/occc_repo/PCA_before_after_batch_correction.png", combined_plot, width = 8, height = 10)
print(combined_plot)

################################################
# rerun limma after batch correction

bio_group <- factor(ifelse(sample_info_sub$Study %in% occc_studies, "OCCC", "ccRCC"))
design <- model.matrix(~ 0 + bio_group)
colnames(design) <- levels(bio_group)

fit <- lmFit(expr_sub, design)

contrast_matrix <- makeContrasts(OCCC_vs_ccRCC = OCCC - ccRCC, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results_corrected <- topTable(fit2, coef = "OCCC_vs_ccRCC", number = Inf, adjust.method = "BH")
#write.csv(results_corrected, file ="/Users/beyzaerkal/Desktop/occc_repo/DE_results_OCCC_vs_ccRCC_batchCorrected.csv", row.names = TRUE)


#######################################
# significant genes for later

fc_thresh <- 1        
fdr_thresh <- 0.05    

results_corrected$Significance <- "Not Sig"
results_corrected$Significance[results_corrected$adj.P.Val < fdr_thresh & results_corrected$logFC > fc_thresh] <- "Up"
results_corrected$Significance[results_corrected$adj.P.Val < fdr_thresh & results_corrected$logFC < -fc_thresh] <- "Down"

sig_genes_corrected <- rownames(results_corrected)[results_corrected$Significance %in% c("Up", "Down")]

results_sig <- results_corrected[sig_genes_corrected, ]
table(results_corrected$Significance)

#write.csv(results_sig, file ="/Users/beyzaerkal/Desktop/occc_repo/DE_genes_OCCC_vs_ccRCC_batchCorrected_significant.csv", row.names = TRUE)
# inside GSEA_2 folder





#batch corr
############ HEATMAP 3
# OC and RC sclaed together


topN <- 500
top_genes <- rownames(results_corrected[order(results_corrected$adj.P.Val), ])[1:topN]
top_genes <- top_genes[top_genes %in% rownames(expr_corrected)]

mat <- expr_corrected[top_genes, , drop = FALSE]
mat_scaled <- t(scale(t(mat)))
mat_scaled[mat_scaled > 3] <- 3
mat_scaled[mat_scaled < -3] <- -3

row_anno <- data.frame(Direction = ifelse(results_corrected[top_genes, "logFC"] > 0, "Up_in_OCCC", "Up_in_ccRCC"))
rownames(row_anno) <- top_genes
row_colors <- c("Up_in_OCCC" = "#e41a1c", "Up_in_ccRCC" = "#377eb8")
ha_row <- rowAnnotation(Direction = row_anno$Direction, col = list(Direction = row_colors), 
                        annotation_name_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                        annotation_legend_param = list(
                          Direction = list(title = "Gene Direction",
                                           title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                                           labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica"))
                        )
)
# OCCC vs ccRCC expression pattern
Heatmap(
  mat_scaled,
  name = "z-score",
  top_annotation = ha_col,
  right_annotation = ha_row,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = function(x) as.dist(1 - cor(t(x), method = "spearman")),
  clustering_distance_columns = function(x) as.dist(1 - cor(x, method = "spearman")),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 7, fontfamily = "Helvetica"),
  row_title = paste0("Top ", topN, " DE genes"), 
  row_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  heatmap_legend_param = list(title = "Expression (Z-score)",
                              title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                              labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")
  )
)

head(results[top_genes, "logFC"], 10)
summary(results[top_genes, "logFC"])


############# HEATMAP 4 (scale separately) - NOT batch corrected (expr_sub)

# balance the top genes to show occc and ccrcc 
topN <- 500
top_up_occc <- head(rownames(results[results$logFC > 0, ]), topN/2)
top_up_ccrcc <- head(rownames(results[results$logFC < 0, ]), topN/2)
top_genes <- c(top_up_occc, top_up_ccrcc)
top_genes <- top_genes[top_genes %in% rownames(expr_sub)] 

mat <- expr_sub[top_genes, , drop = FALSE]

samples_to_plot <- colnames(mat)
rownames(sample_info_sub) <- sample_info_sub$Sample
sample_info_subset <- sample_info_sub[samples_to_plot, , drop = FALSE]
bio_group_subset <- factor(ifelse(sample_info_subset$Study %in% occc_studies, "OCCC", "ccRCC"))
source_type_subset <- factor(sample_info_subset$SourceType)

# remove genes with zero variance in either group - added in here to avoid error 
keep_occc <- apply(mat[, bio_group_subset == "OCCC"], 1, var) > 0
keep_ccrcc <- apply(mat[, bio_group_subset == "ccRCC"], 1, var) > 0
keep_genes <- keep_occc & keep_ccrcc
mat <- mat[keep_genes, ]


# row annot - direction

row_anno <- data.frame(
  Direction = ifelse(results[rownames(mat), "logFC"] > 0, "Up_in_OCCC", "Up_in_ccRCC")
)
rownames(row_anno) <- rownames(mat)
row_anno$Direction <- factor(row_anno$Direction, levels = c("Up_in_OCCC", "Up_in_ccRCC"))
row_colors <- c("Up_in_OCCC" = "#e41a1c", "Up_in_ccRCC" = "#377eb8")
ha_row <- rowAnnotation(Direction = row_anno$Direction, col = list(Direction = row_colors),
                        annotation_name_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                        annotation_legend_param = list(
                          Direction = list(title = "Gene Direction",
                                           title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                                           labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica"))
                        )
)

# scale separtly
mat_scaled <- mat
mat_scaled[, bio_group_subset == "OCCC"] <- t(scale(t(mat[, bio_group_subset == "OCCC"])))
mat_scaled[, bio_group_subset == "ccRCC"] <- t(scale(t(mat[, bio_group_subset == "ccRCC"])))

# cap extremes
mat_scaled[mat_scaled > 3] <- 3
mat_scaled[mat_scaled < -3] <- -3

unique_studies <- unique(sample_info_subset$Study)
study_colors <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(unique_studies)), unique_studies)
bio_colors <- c("OCCC" = "#d73027", "ccRCC" = "#4575b4")
source_colors <- c("cell_line" = "#fdae61", "primary_tumour" = "#313695")

ha_col <- HeatmapAnnotation(
  Study = sample_info_subset$Study,
  TumourType = bio_group_subset,
  SourceType = source_type_subset,
  col = list(
    Study = study_colors,
    TumourType = bio_colors,
    SourceType = source_colors
  ),
  annotation_name_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  annotation_legend_param = list(
    Study = list(title = "Study", ncol = 2, 
                 title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                 labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
    TumourType = list(title = "Tumour Type",
                      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
    SourceType = list(title = "Source Type",
                      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")
    )
  )
)

col_fun <- colorRamp2(c(-3, 0, 3), c("#377eb8", "white", "#e41a1c"))
# OCCC vs ccRCC expression pattern (without Batch-correction)
Heatmap(
  mat_scaled,
  name = "z-score",
  col = col_fun,
  top_annotation = ha_col,
  right_annotation = ha_row,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = function(x) as.dist(1 - cor(t(x), method = "spearman", use = "pairwise.complete.obs")),
  clustering_distance_columns = function(x) as.dist(1 - cor(x, method = "spearman", use = "pairwise.complete.obs")),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 7, fontfamily = "Helvetica"),
  row_title = paste0("Top ", topN, " DE genes"),
  row_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  heatmap_legend_param = list(title = "Expression (Z-score)", 
                              title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                              labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica"))
)


################################
# bacth corrected heatmap 5 (expr_corrected)

# balance the top genes to show occc and ccrcc 
topN <- 500
top_up_occc <- head(rownames(results_corrected[results_corrected$logFC > 0, ]), topN/2)
top_up_ccrcc <- head(rownames(results_corrected[results_corrected$logFC < 0, ]), topN/2)
top_genes <- c(top_up_occc, top_up_ccrcc)
top_genes <- top_genes[top_genes %in% rownames(expr_corrected)] 

mat <- expr_corrected[top_genes, , drop = FALSE]

samples_to_plot <- colnames(mat)
rownames(sample_info_sub) <- sample_info_sub$Sample
sample_info_subset <- sample_info_sub[samples_to_plot, , drop = FALSE]
bio_group_subset <- factor(ifelse(sample_info_subset$Study %in% occc_studies, "OCCC", "ccRCC"))
source_type_subset <- factor(sample_info_subset$SourceType)

# row annot - direction

row_anno <- data.frame(Direction = ifelse(results_corrected[rownames(mat), "logFC"] > 0, "Up_in_OCCC", "Up_in_ccRCC"))
rownames(row_anno) <- rownames(mat)
row_anno$Direction <- factor(row_anno$Direction, levels = c("Up_in_OCCC", "Up_in_ccRCC"))
row_colors <- c("Up_in_OCCC" = "#e41a1c", "Up_in_ccRCC" = "#377eb8")
ha_row <- rowAnnotation(Direction = row_anno$Direction, col = list(Direction = row_colors),
                        annotation_name_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                        annotation_legend_param = list(
                          Direction = list(title = "Gene Direction",
                                           title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                                           labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica"))
                        )
)

# scale separtly
mat_scaled <- mat
mat_scaled[, bio_group_subset == "OCCC"] <- t(scale(t(mat[, bio_group_subset == "OCCC"])))
mat_scaled[, bio_group_subset == "ccRCC"] <- t(scale(t(mat[, bio_group_subset == "ccRCC"])))

# cap extremes
mat_scaled[mat_scaled > 3] <- 3
mat_scaled[mat_scaled < -3] <- -3

unique_studies <- unique(sample_info_subset$Study)
study_colors <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(unique_studies)), unique_studies)
bio_colors <- c("OCCC" = "#d73027", "ccRCC" = "#4575b4")
source_colors <- c("cell_line" = "#fdae61", "primary_tumour" = "#313695")

ha_col <- HeatmapAnnotation(
  Study = sample_info_subset$Study,
  TumourType = bio_group_subset,
  SourceType = source_type_subset,
  col = list(
    Study = study_colors,
    TumourType = bio_colors,
    SourceType = source_colors
  ),
  annotation_legend_param = list(
    Study = list(title = "Study", ncol = 2,
                 title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                 labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
    TumourType = list(title = "Tumour Type",
                      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
    SourceType = list(title = "Source Type",
                      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")
    )
  )
)


col_fun <- colorRamp2(c(-3, 0, 3), c("#377eb8", "white", "#e41a1c"))
# OCCC vs ccRCC expression pattern (Batch-corrected)
Heatmap(
  mat_scaled,
  name = "z-score",
  col = col_fun,
  top_annotation = ha_col,
  right_annotation = ha_row,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = function(x) as.dist(1 - cor(t(x), method = "spearman", use = "pairwise.complete.obs")),
  clustering_distance_columns = function(x) as.dist(1 - cor(x, method = "spearman", use = "pairwise.complete.obs")),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 7, fontfamily = "Helvetica"),
  row_title = paste0("Top ", topN, " DE genes"),
  row_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  heatmap_legend_param = list(title = "Expression (Z-score)",
                              title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                              labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")
  )
)

##################
# HEATMAP 6 batch corrected - without per group scaling
# global scale: no per group scaling to prevent the stripes that is caused after the batch correction

topN <- 500
top_up_occc <- head(rownames(results_corrected[results_corrected$logFC > 0, ]), topN/2)
top_up_ccrcc <- head(rownames(results_corrected[results_corrected$logFC < 0, ]), topN/2)
top_genes <- c(top_up_occc, top_up_ccrcc)
top_genes <- top_genes[top_genes %in% rownames(expr_corrected)]

mat <- expr_corrected[top_genes, , drop = FALSE]

samples_to_plot <- colnames(mat)
sample_info_subset <- sample_info_sub[samples_to_plot, , drop = FALSE]
bio_group_subset <- factor(ifelse(sample_info_subset$Study %in% occc_studies, "OCCC", "ccRCC"))
source_type_subset <- factor(sample_info_subset$SourceType)


row_anno <- data.frame(Direction = ifelse(results_corrected[rownames(mat), "logFC"] > 0, "Up_in_OCCC", "Up_in_ccRCC"))
rownames(row_anno) <- rownames(mat)
row_anno$Direction <- factor(row_anno$Direction, levels = c("Up_in_OCCC", "Up_in_ccRCC"))
row_colors <- c("Up_in_OCCC" = "#e41a1c", "Up_in_ccRCC" = "#377eb8")
ha_row <- rowAnnotation(Direction = row_anno$Direction, col = list(Direction = row_colors),
                        annotation_name_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                        annotation_legend_param = list(
                          Direction = list(title = "Gene Direction",
                                           title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                                           labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica"))
                        )
)


unique_studies <- unique(sample_info_subset$Study)
study_colors <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(unique_studies)), unique_studies)
bio_colors <- c("OCCC" = "#d73027", "ccRCC" = "#4575b4")
source_colors <- c("cell_line" = "#fdae61", "primary_tumour" = "#313695")

ha_col <- HeatmapAnnotation(
  Study = sample_info_subset$Study,
  TumourType = bio_group_subset,
  SourceType = source_type_subset,
  col = list(
    Study = study_colors,
    TumourType = bio_colors,
    SourceType = source_colors
  ),
  annotation_legend_param = list(
    Study = list(title = "Study", ncol = 2, 
                 title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                 labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
    TumourType = list(title = "Tumour Type", 
                      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")),
    SourceType = list(title = "Source Type", 
                      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")
    )
  )
)

# global scaling (center/scale across all samples)
mat_scaled <- t(scale(t(mat), center = TRUE, scale = TRUE))

mat_scaled[mat_scaled > 3] <- 3
mat_scaled[mat_scaled < -3] <- -3

col_fun <- colorRamp2(c(-3, 0, 3), c("#377eb8", "white", "#e41a1c"))
# OCCC vs ccRCC expression pattern (Batch-corrected)
Heatmap(
  mat_scaled,
  name = "z-score",
  col = col_fun,
  top_annotation = ha_col,
  right_annotation = ha_row,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = function(x) as.dist(1 - cor(t(x), method = "spearman", use = "pairwise.complete.obs")),
  clustering_distance_columns = function(x) as.dist(1 - cor(x, method = "spearman", use = "pairwise.complete.obs")),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_names_gp = gpar(fontsize = 7, fontfamily = "Helvetica"),
  row_title = paste0("Top ", topN, " DE genes"),
  row_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  heatmap_legend_param = list(title = "Expression (Z-score)", 
                              title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
                              labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")
  )
)

##################################




##################
# VOLCANO WITHOUT BATCH CORR
##################

fc_thresh <- 1
fdr_thresh <- 0.05

results$Significance <- "Not Sig"
results$Significance[results$adj.P.Val < fdr_thresh & results$logFC > fc_thresh] <- "Up"
results$Significance[results$adj.P.Val < fdr_thresh & results$logFC < -fc_thresh] <- "Down"
results$Significance[results$adj.P.Val == 0 & abs(results$logFC) > 2] <- "Extreme Sig"

# can be added top gene labels on the points to track the extreme sig values (might be technical bias)
top_hits <- results[results$adj.P.Val == 0 & abs(results$logFC) > 2, ]
top_hits$Gene <- rownames(top_hits) # get gene names for labelling

# volcano plot - Volcano Plot: OCCC vs ccRCC
ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_point(data = top_hits, aes(x = logFC, y = -log10(adj.P.Val)), color = "black", size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey", "Extreme Sig" = "black")) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed") +
  geom_text_repel(data = top_hits, aes(x = logFC, y = -log10(adj.P.Val), label = Gene),
                  size = 3, max.overlaps = 20) +
  labs(x = "log2 Fold Change", y = "-log10(FDR)", title = "Volcano Plot: OCCC vs ccRCC") +
  #coord_cartesian(ylim = c(0, 50)) + # cap the y axis by zoom in
  theme_minimal()



common_theme <- theme_minimal(base_size = 7) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    plot.title = element_text(size = 7, face = "bold"),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6)
  )

# volcano plot - black extreme sig data not included
N <- 15
top_up <- head(
  results[results$logFC > 0 & results$adj.P.Val < 0.01, ][order(-results$logFC), ], N)

top_down <- head(
  results[results$logFC < 0 & results$adj.P.Val < 0.01, ][order(results$logFC), ], N)
label_df <- rbind(top_up, top_down)
label_df$Gene <- rownames(label_df)

results <- results[results$Significance != "Extreme Sig", ]

v1 <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed") +
  geom_label_repel(data = label_df, aes(x = logFC, y = -log10(adj.P.Val), label = Gene),
    inherit.aes = FALSE, size = 2.5, max.overlaps = 20, box.padding = 0.35, point.padding = 0.3,
    segment.color = "grey50") +
  labs(x = "log2 Fold Change", y = "-log10(FDR)", title = "Without Batch Correction") +
  common_theme
print(v1)
#########################
# BATCH EFFECT CORR VOLCANO
#########################
results_corrected$Significance <- "Not Sig"
results_corrected$Significance[results_corrected$adj.P.Val < fdr_thresh & results_corrected$logFC > fc_thresh] <- "Up"
results_corrected$Significance[results_corrected$adj.P.Val < fdr_thresh & results_corrected$logFC < -fc_thresh] <- "Down"

table(results_corrected$Significance)

N <- 15
top_up2 <- head(
  results_corrected[results_corrected$logFC > 0 & results_corrected$adj.P.Val < 0.01, ][order(-results_corrected$logFC), ], N)

top_down2 <- head(
  results_corrected[results_corrected$logFC < 0 & results_corrected$adj.P.Val < 0.01, ][order(results_corrected$logFC), ], N)
label_df2 <- rbind(top_up2, top_down2)
label_df2$Gene <- rownames(label_df2)


v2 <- ggplot(results_corrected, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed") +
  geom_label_repel(data = label_df2, aes(x = logFC, y = -log10(adj.P.Val), label = Gene),
                   inherit.aes = FALSE, size = 2.5, max.overlaps = 20, box.padding = 0.35, point.padding = 0.3,
                   segment.color = "grey50") +
  labs(x = "log2 Fold Change", y = "-log10(FDR)", title = "With Batch Correction") +
  common_theme
print(v2)
# combine both
combined_volcano <- plot_grid(v1, v2, ncol = 2)
print(combined_volcano)
#ggsave("/Users/beyzaerkal/Desktop/occc_repo/Volcano_plots_before_after_batch_correction.png", combined_volcano, width = 12, height = 6)



####################
# SAVE COMBINED FIGURES
####################
# PCA + Volcano + Heatmap

volcano_block <- wrap_plots(v1, v2, ncol = 2)
heatmap_block <- wrap_elements(full = as.ggplot(h2))

layout_plot <- wrap_plots(
  volcano_block,
  heatmap_block,
  ncol = 1,
  heights = c(0.75, 1.3)
)

rna_vol_heatmap_plot <- ggdraw(layout_plot) +
  draw_plot_label(
    label = c("A", "B"),
    x = 0.01,
    y = c(0.97, 0.50),
    hjust = 0,
    vjust = 1,
    size = 7,
    fontface = "bold"
  )

print(rna_vol_heatmap_plot)
# 7 x 8

####################
# input the rds prot heatmap and combine with panel purposes
hm_inputs <- readRDS("/Users/beyzaerkal/Desktop/occc_repo/DEProt_heatmap_inputs.RDS")

mat_prot <- hm_inputs$mat
sig_genes_prot <- hm_inputs$sig_genes
sample_type <- hm_inputs$sample_type

row_type <- ifelse(sig_genes_prot[rownames(mat_prot), "logFC"] > 0,
                   "Up in OCCC", "Up in ccRCC")
row_type <- factor(row_type, levels = c("Up in OCCC", "Up in ccRCC"))

ha_row <- rowAnnotation(
  Direction = row_type,
  col = list(Direction = c("Up in OCCC" = "#e41a1c",
                           "Up in ccRCC" = "#377eb8")),
  annotation_name_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  annotation_legend_param = list(
    Direction = list(
      title = "Gene Direction",
      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")
    )
  )
)

ha_col <- HeatmapAnnotation(
  TumourType = sample_type,
  col = list(TumourType = c("OCCC" = "#e41a1c",
                            "ccRCC" = "#377eb8")),
  annotation_name_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  annotation_legend_param = list(
    TumourType = list(
      title = "Tumour Type",
      title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
      labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")
    )
  )
)

col_fun <- colorRamp2(c(-3, 0, 3), c("#377eb8", "white", "#e41a1c"))

h_prot <- Heatmap(
  mat_prot,
  name = "z-score",
  col = col_fun,
  left_annotation = ha_row,
  top_annotation = ha_col,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_title = "Top 100 DE genes",
  row_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
  heatmap_legend_param = list(
    title = "Expression (Z-score)",
    title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
    labels_gp = gpar(fontsize = 6, fontfamily = "Helvetica")
  )
)
print(h_prot)

# COMBINED PANELS

heatmap_block_prot <- draw(h_prot)
h_prot_grob <- grid.grabExpr(draw(h_prot))
volcano_block <- wrap_plots(v1, v2, ncol = 2)
heatmap_block_rna <- wrap_elements(full = as.ggplot(h2))

layout_plot <- wrap_plots(
  volcano_block,
  heatmap_block_rna,
  h_prot_grob,
  ncol = 1,
  heights = c(0.75, 1.3, 1)
)

rna_vol_heatmap_plot <- ggdraw(layout_plot) +
  draw_plot_label(
    label = c("A", "B", "C"),
    x = 0.01,
    y = c(0.97, 0.65, 0.30),
    hjust = 0,
    vjust = 1,
    size = 7,
    fontface = "bold"
  )

print(rna_vol_heatmap_plot)


# COMBINED PANELS 2 wihto only batch corrected

heatmap_block_prot <- draw(h_prot)
h_prot_grob <- grid.grabExpr(draw(h_prot))
volcano_pca_block <- wrap_plots(pca_after_plot, v2, ncol = 2)
heatmap_block_rna <- wrap_elements(full = as.ggplot(h2))

layout_plot <- wrap_plots(
  volcano_pca_block,
  heatmap_block_rna,
  h_prot_grob,
  ncol = 1,
  heights = c(0.75, 1.3, 1)
)

rna_vol_heatmap_plot <- ggdraw(layout_plot) +
  draw_plot_label(
    label = c("A", "B", "C", "D"),
    x = c(0.01, 0.51, 0.01, 0.01),
    y = c(0.97, 0.97, 0.65, 0.30),
    hjust = 0,
    vjust = 1,
    size = 7,
    fontface = "bold"
  )

print(rna_vol_heatmap_plot)


##################
# ORA - clusterprofiler ""without"" batch correction

sig_genes <- results[results$adj.P.Val < 0.05, ]
gene_symbols <- rownames(sig_genes)
# cnv to entrez
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) 
# merge w DE
sig_genes$SYMBOL <- rownames(sig_genes)
merged <- merge(entrez_ids, sig_genes, by = "SYMBOL")


# GO enrich - bio process
go_bp <- enrichGO(gene = merged$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",  
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# Go cellular components
go_cc <- enrichGO(gene = merged$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",      
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# GO molecular function
go_mf <- enrichGO(gene = merged$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF", 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)


# kegg
kegg <- enrichKEGG(gene = merged$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)


dotplot(go_bp, showCategory = 20, title = "GO Biological Processes")
dotplot(go_cc, showCategory = 20) + ggtitle("GO Cellular Component")
dotplot(go_mf, showCategory = 20) + ggtitle("GO Molecular Function")
barplot(kegg, showCategory = 20, color = "p.adjust") + ggtitle("KEGG Pathway Enrichment") + theme_minimal()

write.csv(as.data.frame(go_bp), "GO_BP_enrichment.csv")
write.csv(as.data.frame(kegg), "KEGG_enrichment.csv")


reactome <- enrichPathway(gene = merged$ENTREZID, organism = "human", pvalueCutoff = 0.05)
write.csv(as.data.frame(reactome), "REACTOME_enrichment.csv")

# multiple enrichment analysis #
# combine all

go_bp_df <- as.data.frame(go_bp)
go_bp_df$Source <- "GO_BP"

kegg_df <- as.data.frame(kegg)
kegg_df$Source <- "KEGG"

reactome_df <- as.data.frame(reactome)
reactome_df$Source <- "Reactome"


common_cols <- Reduce(intersect, list(colnames(go_bp_df), colnames(kegg_df), colnames(reactome_df)))

go_bp_df <- go_bp_df[, c(common_cols, "Source")]
kegg_df <- kegg_df[, c(common_cols, "Source")]
reactome_df <- reactome_df[, c(common_cols, "Source")]


combined <- rbind(go_bp_df, kegg_df, reactome_df)


combined <- combined %>%
  group_by(Source) %>%
  slice_min(order_by = p.adjust, n = 20) %>%
  ungroup()

combined$Description <- factor(combined$Description, levels = rev(unique(combined$Description)))

ggplot(combined, aes(x = Source, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Multi-Library Functional Enrichment",
       x = "Gene Set Library",
       y = "Enriched Term",
       size = "Gene Count",
       color = "Adjusted P-value") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8))


#######################
# Batch Corrected
#######################
# ranked GSEA
# using t stats

ranked_genes <- topTable(fit2, number = Inf)
gene_symbols <- rownames(ranked_genes)

symbol2entrez <- bitr(gene_symbols,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)


mapped_genes <- ranked_genes[rownames(ranked_genes) %in% symbol2entrez$SYMBOL, ]
entrez_ids <- symbol2entrez$ENTREZID[match(rownames(mapped_genes), symbol2entrez$SYMBOL)]

# ranked vector
gene_list <- mapped_genes$t   # logFC/t
names(gene_list) <- entrez_ids
gene_list <- gene_list[!is.na(names(gene_list))] # unmapped removed
gene_list <- sort(gene_list, decreasing = TRUE)
gsea_go <- gseGO(geneList = gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.2,
                 verbose = FALSE)

gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = "hsa",
                     minGSSize = 10,
                     pvalueCutoff = 0.2) 

gsea_reactome <- gsePathway(geneList = gene_list,
                            organism = "human",
                            minGSSize = 10,
                            pvalueCutoff = 0.05)

dotplot(gsea_go, showCategory = 20, color = "NES") + scale_color_gradient(low = "blue", high = "red") + ggtitle("GSEA GO Biological Process") 
dotplot(gsea_kegg, showCategory = 20, color = "NES") + scale_color_gradient(low = "blue", high = "red") + ggtitle("GSEA KEGG Pathways")
dotplot(gsea_reactome, showCategory = 20, color = "NES") + scale_color_gradient(low = "blue", high = "red") + ggtitle("GSEA Reactome Pathways")

gsea_kegg_net <- pairwise_termsim(gsea_kegg)
emapplot(gsea_kegg_net)

gseaplot2(gsea_reactome, geneSetID = 1, title = gsea_reactome$Description[1], color = "#7570b3") 

# select the relevant gsea_reactome plots
gseaplot2(gsea_reactome, geneSetID = "R-HSA-2428928", title = "IRS-related events triggered by IGF1R", color = "#7570b3")
# main figure - added to the panel in GSEA_NES_heatmap file
gseaplot2(gsea_reactome, geneSetID = "R-HSA-2404192", title = "IGF1R signaling cascade", color = "#1b9e77")
gseaplot2(gsea_reactome, geneSetID = "R-HSA-912694", title = "Regulation of IFNA/IFNB signaling", color = "#d95f02")
gseaplot2(gsea_reactome, geneSetID = "R-HSA-5654743", title = "Signaling by FGFR4", color = "#e7298a")


#require(DOSE)
dotplot(gsea_reactome, showCategory = 10, split = ".sign") + facet_grid(.~.sign) + theme(text = element_text(family = "Helvetica")) +
  theme(legend.title = element_text(size = 7),  
        axis.text.x = element_text(size = 7),  
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 7, face = "bold"))

#write.csv(as.data.frame(gsea_go), file = "/Users/beyzaerkal/Desktop/occc_repo/RNA_GSEA_GO_BP.csv")
#write.csv(as.data.frame(gsea_kegg), file = "/Users/beyzaerkal/Desktop/occc_repo/RNA_GSEA_KEGG.csv")
#write.csv(as.data.frame(gsea_reactome), file = "/Users/beyzaerkal/Desktop/occc_repo/RNA_GSEA_Reactome.csv")
saveRDS(gsea_reactome, file = "/Users/beyzaerkal/Desktop/occc_repo/RNA_GSEA_Reactome.rds")


# summary - top pathways
head(gsea_go@result)

# filtered gsea
exclude_keywords <- c("sensory", "olfactory", "smell", "vision", "taste", "visual", "stimulus", "sperm", "keratin")
head(gsea_reactome@result$Description)
gsea_reactome_filtered <- gsea_reactome
gsea_reactome_filtered@result <- gsea_reactome@result[!grepl(paste(exclude_keywords, collapse="|"),
                                                       gsea_reactome@result$Description, ignore.case = TRUE), ]

gsea_reactome_filtered@result <- gsea_reactome_filtered@result[order(gsea_reactome_filtered@result$p.adjust), ]

head(gsea_reactome_filtered, 10) # by NES
dotplot(gsea_reactome_filtered, showCategory = 20, color = "NES") +
  scale_color_gradient(low = "blue", high = "red") +
  ggtitle("Filtered GSEA Reactome Pathways")

dotplot(gsea_reactome_filtered, showCategory = 10, split = ".sign") + facet_grid(.~.sign) # by p-adj significance

# saved as: Filtered Reactome GSEA dotplot
dotplot(gsea_reactome_filtered, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) + theme(text = element_text(family = "Helvetica")) +
  theme(legend.title = element_text(size = 7),  
        axis.text.x = element_text(size = 7),  
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        strip.text = element_text(size = 7, face = "bold"))


########################
# KINASES AND ONCOGENES
########################


#  oncogene list 
oncokb <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/cancerGeneList.tsv", show_col_types = FALSE)

# clean col
oncokb_clean <- oncokb[, c("Hugo Symbol", "Gene Type", "OncoKB Annotated")]
colnames(oncokb_clean)[1] <- "Gene"  # rename 

# adding gene names to the sig_genes dataframe as a column
sig_genes$Gene <- rownames(sig_genes)
# merge DE results 
sig_annotated <- merge(sig_genes, oncokb_clean, by = "Gene", all.x = TRUE)
head(sig_annotated)
# Only keep DE genes that are oncogenes or tumor suppressors
#sig_onco_tsg <- sig_annotated[!is.na(sig_annotated$`Gene Type`), ]
sig_onco_tsg <- sig_annotated[sig_annotated$`Gene Type` %in% c("ONCOGENE", "TSG", "ONCOGENE_AND_TSG"), ]

#write.csv(sig_onco_tsg, "significant_oncogenes_tsgs.csv", row.names = FALSE)

ggplot(sig_onco_tsg, aes(x = logFC, y = -log10(adj.P.Val), color = `Gene Type`)) +
  geom_point() +
  geom_text(aes(label = Gene), vjust = 1, hjust = 1, size = 3, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = "Differential Expression of Oncogenes and Tumour Suppressor Genes", x = "log2 Fold Change", y = "-log10(FDR)")



# significant oncogenes clusterprofiler
oncogenes <- oncokb %>% filter(`Gene Type` == "ONCOGENE") %>% pull(`Hugo Symbol`) %>% toupper()
# to get sig_genes_filtered run kinases first
# reorder and filter significant genes with logFC threshold
sig_genes_filtered <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
sig_genes_filtered$Gene <- toupper(rownames(sig_genes_filtered))

sig_oncogenes <- sig_genes_filtered %>% filter(Gene %in% oncogenes)
# Oncogenes
oncogene_entrez <- bitr(sig_oncogenes$Gene, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Oncogenes
ego_onco <- enrichGO(gene = oncogene_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, readable = TRUE)

# Oncogenes kegg
ekegg_onco <- enrichKEGG(gene = oncogene_entrez$ENTREZID,
                         organism = 'hsa', pvalueCutoff = 0.05)


reactome_onco <- enrichPathway(gene = oncogene_entrez$ENTREZID, organism = "human", pvalueCutoff = 0.05)

# separate plots
barplot(ego_onco, showCategory = 10, title = "GO Enrichment - Oncogenes")
dotplot(ego_onco, showCategory = 10, title = "GO Dotplot - Oncogenes")
dotplot(ekegg_onco, showCategory = 10, title = "KEGG - Oncogenes")
dotplot(reactome_onco, showCategory = 10, title = "REACTOME - Oncogenes")

cnetplot(ego_onco, categorySize = "pvalue")


# combined oncogene plot ( threshold 0.01)

go_bp_onco_df <- as.data.frame(ego_onco)
go_bp_onco_df$Source <- "GO_BP"

kegg_onco_df <- as.data.frame(ekegg_onco)
kegg_onco_df$Source <- "KEGG"

reactome_onco_df <- as.data.frame(reactome_onco)
reactome_onco_df$Source <- "Reactome"


common_cols <- Reduce(intersect, list(colnames(go_bp_onco_df), colnames(kegg_onco_df), colnames(reactome_onco_df)))

go_bp_onco_df <- go_bp_onco_df[, c(common_cols, "Source")]
kegg_onco_df <- kegg_onco_df[, c(common_cols, "Source")]
reactome_onco_df <- reactome_onco_df[, c(common_cols, "Source")]

combined_onco <- rbind(go_bp_onco_df, kegg_onco_df, reactome_onco_df)
# 0.01 too many pathways
filtered_onco <- combined_onco %>% filter(p.adjust < 0.01)


ggplot(filtered_onco, aes(x = Source, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = "Multi-Library Functional Enrichment ONCOGENE",
       x = "Gene Set Library",
       y = "Enriched Term",
       size = "Gene Count",
       color = "Adjusted P-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


##########################
# kinase list ( threshold 0.05)
kinases <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/kinase_basic.xlsx", col_names = TRUE)
colnames(kinases)

kinases$Offical_gene_symbol <- toupper(kinases$Offical_gene_symbol)

sig_genes_k <- results_corrected[results_corrected$adj.P.Val < 0.05, ]
sig_genes_k$Gene <- toupper(rownames(sig_genes_k))

# significant genes -- kinase list
sig_kinases <- sig_genes_k[sig_genes_k$Gene %in% kinases$Offical_gene_symbol, ]

head(sig_kinases) # 401   7

#write.csv(sig_genes_k, "significant_kinases_RNA.csv", row.names = FALSE)
###

kinase_entrez <- bitr(sig_kinases$Gene,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

go_kinase <- enrichGO(gene = kinase_entrez$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",  # other "MF" "CC" 
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

kegg_kinase <- enrichKEGG(gene = kinase_entrez$ENTREZID,
                           organism = 'hsa',
                           pvalueCutoff = 0.05)

reactome_kinase <- enrichPathway(gene = kinase_entrez$ENTREZID,
                                  organism = "human",
                                  pvalueCutoff = 0.05,
                                  readable = TRUE)



dotplot(go_kinase, showCategory = 10, title = "GO Enrichment - Kinases")
dotplot(kegg_kinase, showCategory = 10, title = "KEGG Pathways - Kinases")
dotplot(reactome_kinase, showCategory = 20, title = "Reactome Pathways - Kinases")


cnetplot(reactome_kinase, categorySize = "pvalue")

# combined kinase plot ( threshold 0.05)

go_kinase_df <- as.data.frame(go_kinase)
go_kinase_df$Source <- "GO_BP"

kegg_kinase_df <- as.data.frame(kegg_kinase)
kegg_kinase_df$Source <- "KEGG"

reactome_kinase_df <- as.data.frame(reactome_kinase)
reactome_kinase_df$Source <- "Reactome"


common_cols <- Reduce(intersect, list(colnames(go_kinase_df), colnames(kegg_kinase_df), colnames(reactome_kinase_df)))

go_kinase_df <- go_kinase_df[, c(common_cols, "Source")]
kegg_kinase_df <- kegg_kinase_df[, c(common_cols, "Source")]
reactome_kinase_df <- reactome_kinase_df[, c(common_cols, "Source")]

combined_kinase <- rbind(go_kinase_df, kegg_kinase_df, reactome_kinase_df)
# 0.01 if too many pathways
filtered_kinase <- combined_kinase %>% filter(p.adjust < 0.01)


ggplot(filtered_kinase, aes(x = Source, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = "Multi-Library Functional Enrichment KINASES",
       x = "Gene Set Library",
       y = "Enriched Term",
       size = "Gene Count",
       color = "Adjusted P-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


#library(ggnewscale)


##################################################################
# Batch corrected ORA (2)
##################################################################

# ORA - clusterprofiler

sig_genes <- results_corrected[results_corrected$adj.P.Val < 0.05, ]
gene_symbols <- rownames(sig_genes)

entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) 

# merge w DE
sig_genes$SYMBOL <- rownames(sig_genes)
merged <- merge(entrez_ids, sig_genes, by = "SYMBOL")


# GO enrich - bio process
go_bp <- enrichGO(gene = merged$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",  
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# Go cellular components
go_cc <- enrichGO(gene = merged$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",      
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# GO molecular function
go_mf <- enrichGO(gene = merged$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF", 
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)


# kegg
kegg <- enrichKEGG(gene = merged$ENTREZID,
                   organism = 'hsa',
                   pvalueCutoff = 0.05)


dotplot(go_bp, showCategory = 20, title = "GO Biological Process")
dotplot(go_cc, showCategory = 20) + ggtitle("GO Cellular Component")
dotplot(go_mf, showCategory = 20) + ggtitle("GO Molecular Function")
barplot(kegg, showCategory = 20, color = "p.adjust") + ggtitle("KEGG Pathway Enrichment") + theme_minimal()

write.csv(as.data.frame(go_bp), "GO_BP_enrichment.csv")
write.csv(as.data.frame(kegg), "KEGG_enrichment.csv")


reactome <- enrichPathway(gene = merged$ENTREZID, organism = "human", pvalueCutoff = 0.05)

write.csv(as.data.frame(reactome), "REACTOME_enrichment.csv")

# multiple enrichment analysis #
# combine all

go_bp_df <- as.data.frame(go_bp)
go_bp_df$Source <- "GO_BP"

kegg_df <- as.data.frame(kegg)
kegg_df$Source <- "KEGG"

reactome_df <- as.data.frame(reactome)
reactome_df$Source <- "Reactome"


common_cols <- Reduce(intersect, list(colnames(go_bp_df), colnames(kegg_df), colnames(reactome_df)))

go_bp_df <- go_bp_df[, c(common_cols, "Source")]
kegg_df <- kegg_df[, c(common_cols, "Source")]
reactome_df <- reactome_df[, c(common_cols, "Source")]


combined <- rbind(go_bp_df, kegg_df, reactome_df)


combined <- combined %>%
  group_by(Source) %>%
  slice_min(order_by = p.adjust, n = 20) %>%
  ungroup()

combined$Description <- factor(combined$Description, levels = rev(unique(combined$Description)))

ggplot(combined, aes(x = Source, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Multi-Library Functional Enrichment",
       x = "Gene Set Library",
       y = "Enriched Term",
       size = "Gene Count",
       color = "Adjusted P-value") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8))



###############################################
# batch corrected KINASES and ONCOGENES
###############################################

# after batch correction

sig_genes <- read_csv("/Users/beyzaerkal/Desktop/occc_repo/DE_genes_OCCC_vs_ccRCC_batchCorrected_significant.csv")
colnames(sig_genes)[1] <- "Gene"

# oncogene list 
oncokb <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/cancerGeneList.tsv", show_col_types = FALSE)
oncokb_clean <- oncokb[, c("Hugo Symbol", "Gene Type", "OncoKB Annotated")]
colnames(oncokb_clean)[1] <- "Gene" 

# merge DE results with OncoKB
sig_annotated <- merge(sig_genes, oncokb_clean, by = "Gene", all.x = TRUE)
head(sig_annotated)

# filter for oncogenes and TSGs
sig_onco_tsg <- sig_annotated %>% filter(`Gene Type` %in% c("ONCOGENE", "TSG", "ONCOGENE_AND_TSG"))


write.csv(sig_onco_tsg, "/Users/beyzaerkal/Desktop/occc_repo/oncogenes/significant_oncogenes_tsgs.csv", row.names = FALSE)

ggplot(sig_onco_tsg, aes(x = logFC, y = -log10(adj.P.Val), color = `Gene Type`)) +
  geom_point() +
  geom_text(aes(label = Gene), vjust = 1, hjust = 1, size = 3, check_overlap = TRUE) +
  theme_minimal() +
  labs(title = "Differential Expression of Oncogenes and Tumour Suppressor Genes", x = "log2 Fold Change", y = "-log10(FDR)")

ggplot(sig_onco_tsg, aes(x = logFC, y = -log10(adj.P.Val), color = `Gene Type`)) +
  geom_point() +
  geom_text(aes(label = Gene), vjust = 1, hjust = 1, size = 3, check_overlap = TRUE) +
geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +          
  geom_vline(xintercept = -1, linetype = "dashed", color = "black") +         
  theme_minimal() +
  labs(title = "Differential Expression of Oncogenes and Tumour Suppressor Genes",
    x = "log2 Fold Change",
    y = "-log10(FDR)")

# significant oncogenes only 
sig_oncogenes <- sig_onco_tsg %>% filter(`Gene Type` == "ONCOGENE")

# Oncogenes
oncogene_entrez <- bitr(sig_oncogenes$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Oncogenes
ego_onco <- enrichGO(gene = oncogene_entrez$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, readable = TRUE)

# Oncogenes kegg
ekegg_onco <- enrichKEGG(gene = oncogene_entrez$ENTREZID,
                         organism = 'hsa', pvalueCutoff = 0.05)


reactome_onco <- enrichPathway(gene = oncogene_entrez$ENTREZID, organism = "human", pvalueCutoff = 0.05)

# separate plots
dotplot(ego_onco, showCategory = 10, title = "GO BP - Oncogenes")
barplot(ekegg_onco, showCategory = 10, title = "KEGG - Oncogenes")
barplot(reactome_onco, showCategory = 10, title = "REACTOME - Oncogenes")
barplot(reactome_onco, showCategory = 20, title = "REACTOME - Oncogenes")

cnetplot(ego_onco, categorySize = "pvalue", circular = TRUE, color_category = "black")


# combined oncogene plot ( threshold 0.01)

go_bp_onco_df <- as.data.frame(ego_onco)
go_bp_onco_df$Source <- "GO_BP"

kegg_onco_df <- as.data.frame(ekegg_onco)
kegg_onco_df$Source <- "KEGG"

reactome_onco_df <- as.data.frame(reactome_onco)
reactome_onco_df$Source <- "Reactome"


common_cols <- Reduce(intersect, list(colnames(go_bp_onco_df), colnames(kegg_onco_df), colnames(reactome_onco_df)))

go_bp_onco_df <- go_bp_onco_df[, c(common_cols, "Source")]
kegg_onco_df <- kegg_onco_df[, c(common_cols, "Source")]
reactome_onco_df <- reactome_onco_df[, c(common_cols, "Source")]

combined_onco <- rbind(go_bp_onco_df, kegg_onco_df, reactome_onco_df)

top_n <- 30 # also 20 applied
filtered_onco <- combined_onco %>% group_by(Source) %>% slice_min(order_by = p.adjust, n = top_n) %>% ungroup()

ggplot(filtered_onco, aes(x = Source, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = "Multi-Library Functional Enrichment ONCOGENE",
       x = "Gene Set Library",
       y = "Enriched Term",
       size = "Gene Count",
       color = "Adjusted P-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


##########################
# kinase list ( threshold 0.05)
kinases <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/kinase_basic.xlsx", col_names = TRUE)
colnames(kinases)

kinases$Offical_gene_symbol <- toupper(kinases$Offical_gene_symbol)

# significant genes -- kinase list
sig_kinases <- sig_genes %>% dplyr::filter(Gene %in% kinases$Offical_gene_symbol)

head(sig_kinases) 
nrow(sig_kinases)

#write.csv(sig_kinases, "/Users/beyzaerkal/Desktop/occc_repo/kinases/significant_kinases.csv", row.names = FALSE)
###

kinase_entrez <- bitr(sig_kinases$Gene,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)

go_kinase <- enrichGO(gene = kinase_entrez$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      ont = "BP",  # other "MF" "CC" 
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      readable = TRUE)

kegg_kinase <- enrichKEGG(gene = kinase_entrez$ENTREZID,
                          organism = 'hsa',
                          pvalueCutoff = 0.05)

reactome_kinase <- enrichPathway(gene = kinase_entrez$ENTREZID,
                                 organism = "human",
                                 pvalueCutoff = 0.05,
                                 readable = TRUE)


barplot(go_kinase, showCategory = 10, title = "GO BP - Kinases")
barplot(kegg_kinase, showCategory = 10, title = "KEGG Pathways - Kinases")
barplot(reactome_kinase, showCategory = 10, title = "Reactome Pathways - Kinases")

# to combine images in anotherh script
saveRDS(kegg_kinase,      "rna_kegg_kinase.rds")
saveRDS(reactome_kinase,  "rna_reactome_kinase.rds")

cnetplot(reactome_kinase, categorySize = "pvalue")

# combined kinase plot ( threshold 0.05)

go_kinase_df <- as.data.frame(go_kinase)
go_kinase_df$Source <- "GO_BP"

kegg_kinase_df <- as.data.frame(kegg_kinase)
kegg_kinase_df$Source <- "KEGG"

reactome_kinase_df <- as.data.frame(reactome_kinase)
reactome_kinase_df$Source <- "Reactome"


common_cols <- Reduce(intersect, list(colnames(go_kinase_df), colnames(kegg_kinase_df), colnames(reactome_kinase_df)))

go_kinase_df <- go_kinase_df[, c(common_cols, "Source")]
kegg_kinase_df <- kegg_kinase_df[, c(common_cols, "Source")]
reactome_kinase_df <- reactome_kinase_df[, c(common_cols, "Source")]

combined_kinase <- rbind(go_kinase_df, kegg_kinase_df, reactome_kinase_df)
# 0.01 if too many pathways
filtered_kinase <- combined_kinase %>% filter(p.adjust < 0.01)

top_n <- 30 # also 20 applied
filtered_combined_kinase <- combined_kinase %>% group_by(Source) %>% slice_min(order_by = p.adjust, n = top_n) %>% ungroup()

ggplot(filtered_combined_kinase, aes(x = Source, y = Description, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.7) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title = "Multi-Library Functional Enrichment KINASES",
       x = "Gene Set Library",
       y = "Enriched Term",
       size = "Gene Count",
       color = "Adjusted P-value") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

#####################
# combine kinase reactome + prot reactome
#####################


prot_kegg_kinase <- readRDS("prot_kegg_kinase.rds")
prot_reactome_kinase <- readRDS("prot_reactome_kinase.rds")

p_kegg_prot_kinase <- dotplot(prot_kegg_kinase, showCategory = 10)
p_reactome_prot_kinase <- dotplot(prot_reactome_kinase, showCategory = 10)

print(p_kegg_prot_kinase)
print(p_reactome_prot_kinase)

# rna
p_kegg_rna_kinase <- dotplot(kegg_kinase, showCategory = 10, title = "KEGG Pathways - Kinases")
p_reactome_rna_kinase <- dotplot(reactome_kinase, showCategory = 10, title = "Reactome Pathways - Kinases")


# main
p_reactome_prot_kinase | p_reactome_rna_kinase

#supp material
p_kegg_prot_kinase | p_kegg_rna_kinase



#############
# RNA MSigDB C3 ONLY
#############

DE_RNA <- read.csv("/Users/beyzaerkal/Desktop/occc_repo/DE_results_OCCC_vs_ccRCC_batchCorrected.csv", check.names = FALSE)
gene_sets_path <- "/Users/beyzaerkal/Desktop/internship/internship_env/c3.all.v2025.1.Hs.symbols.gmt"
colnames(DE_RNA)[1] <- "Geneid"

gene_sets_TFT <- read.gmt(gene_sets_path)

unique_terms <- unique(gene_sets_TFT$term)
length(unique_terms)

# convert all to uppercase
gene_sets_TFT$gene <- toupper(gene_sets_TFT$gene)

de_genes <- DE_RNA[abs(DE_RNA$logFC) > 0.5, ]


prot_de_symbols <- toupper(de_genes$Geneid)
universe_genes  <- toupper(DE_RNA$Geneid)


#custome ORA
ora_C3 <- enricher(gene = prot_de_symbols,
                   universe = universe_genes,
                   TERM2GENE = gene_sets_TFT,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

go_C3MF <- enrichGO(gene = prot_de_symbols,
                    OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05)

sum(prot_de_symbols %in% gene_sets_TFT$gene)
dotplot(ora_C3, showCategory = 15, title = "")
dotplot(go_C3MF, showCategory = 15, title = "GO MF - C3 Genes")
# stats sig 

ora_df <- as.data.frame(ora_C3)
sig_C3terms <- subset(ora_df, p.adjust < 0.05)

head(sig_C3terms)
write.csv(sig_C3terms, "RNA_significant_C3_terms.csv", row.names = FALSE)

