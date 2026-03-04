
#################
# RNA-seq only
#################

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(limma)
library(readxl)
library(biomaRt)
library(openxlsx)


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
# expression matrix of both ovary and renal
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

# remove extra col info from gtex
ovary_gtex_filtered <- ovary_gtex_filtered[, !colnames(ovary_gtex_filtered) %in% "gene_biotype"]
renal_gtex_filtered <- renal_gtex_filtered[, !colnames(renal_gtex_filtered) %in% "gene_biotype"]

# combine RAW TPM first (NOT log)
common_genes_all <- Reduce(intersect, list(
  rownames(occc_all),
  rownames(ovary_gtex_filtered),
  rownames(tpm_ccRCC),
  rownames(renal_gtex_filtered)
))

length(common_genes_all)

# group all datasets to common genes
occc_all <- occc_all[common_genes_all, ]
ovary_gtex_filtered <- ovary_gtex_filtered[common_genes_all, ]
tpm_ccRCC <- tpm_ccRCC[common_genes_all, ]
renal_gtex_filtered <- renal_gtex_filtered[common_genes_all, ]

# combine raw tpm
expr_raw_oc_rc <- cbind(occc_all, ovary_gtex_filtered, tpm_ccRCC, renal_gtex_filtered)
dim(expr_raw_oc_rc)

# filter low expression genes
# keep only genes with tpm >= 1 in at least 6 samples
keep_genes <- rowSums(expr_raw_oc_rc >= 1) >= 6
expr_raw_oc_rc_filtered <- expr_raw_oc_rc[keep_genes, ]
nrow(expr_raw_oc_rc) # 15488
nrow(expr_raw_oc_rc_filtered) # All 15,488 common genes already have tpm >= 1 in at least 6 samples 

# log2 transform 
# careful naming on from expr_raw to expr_mat for DE analysis
expr_mat_oc_rc_filtered <- log2(expr_raw_oc_rc_filtered + 1)

# filter low variance genes - remove genes with zero variance across all samples
expr_mat_oc_rc_filtered <- expr_mat_oc_rc_filtered[apply(expr_mat_oc_rc_filtered, 1, var) > 0, ]
nrow(expr_mat_oc_rc_filtered) # 15488

#saveRDS(expr_mat_oc_rc_filtered, file = "/Users/beyzaerkal/Desktop/occc_multi-omics/processed/expression_matrix_occc_renal.rds")

##############
samples <- colnames(expr_mat_oc_rc_filtered)

# GTEx ovary and renal sample IDs - for separate labelling

ovary_samples <- colnames(ovary_gtex_filtered)
renal_samples <- colnames(renal_gtex_filtered)
ccrcc_samples <- colnames(tpm_ccRCC)


ea_samples <- c("ES-2", "JHOC-5", "OVISE", "OVMANA", "OVTOKO", "RMG-I", "TOV-21G")

study <- ifelse(samples %in% ea_samples, "Expression Atlas (2026)",
                ifelse(samples %in% c("C1", "C2", "C3", "C4", "C5", "C6"), "Nagasawa et al. (2019)",
                       ifelse(grepl("^CCC_", samples), "GSE189553",
                              ifelse(grepl("^OVA[0-9]+", samples), "GSE160692",
                                     ifelse(grepl("IGO_07456", samples), "Bolton et al. (2022)",
                                            ifelse(samples %in% ccrcc_samples, "TCGA-KIRC",
                                                   ifelse(samples %in% ovary_samples, "GTEx_Ovary",
                                                          ifelse(samples %in% renal_samples, "GTEx_Renal_Cortex", NA))))))))

# map tisseu form study
tissue <- ifelse(study %in% c("Expression Atlas (2026)", "Nagasawa et al. (2019)",
                              "GSE189553", "GSE160692", "Bolton et al. (2022)", "GTEx_Ovary"), "Ovary",
                 ifelse(study %in% c("TCGA-KIRC", "GTEx_Renal_Cortex"), "Kidney", NA))

status <- ifelse(samples %in% ovary_samples | samples %in% renal_samples, "Normal", "Tumour")
source_type <- ifelse(study == "Expression Atlas (2026)", "Cell line", "Primary Tumour")

# metadata combined
sample_info <- data.frame(
  Sample = samples,
  Study = factor(study),
  Tissue = factor(tissue),
  Status = factor(status),
  SourceType = factor(source_type)
)

#saveRDS(sample_info, file = "/Users/beyzaerkal/Desktop/occc_multi-omics/processed/sample_info_OCandRC.rds")

rownames(sample_info) <- sample_info$Sample
samples[duplicated(samples)] # col checks that have duplicates
all(colnames(expr_mat_oc_rc_filtered) == rownames(sample_info)) 

# Combined group factor - make categorical data combining metadata into groups
group <- rep(NA, length(samples))
names(group) <- samples
group[samples %in% colnames(occc_all)] <- "OCCC"
group[samples %in% colnames(ovary_gtex_filtered)] <- "GTEx_Ovary"
group[samples %in% colnames(tpm_ccRCC)] <- "TCGA-KIRC" # ccRCC
group[samples %in% colnames(renal_gtex_filtered)] <- "GTEx_Renal_Cortex"

sample_info$Group <- factor(group, levels = c("GTEx_Ovary", "OCCC", "GTEx_Renal_Cortex", "TCGA-KIRC"))

table(sample_info$Group)
table(sample_info$Tissue, sample_info$Status)

##############################
# gene checks
# Use top 500 most variable genes 
top_genes <- head(order(apply(expr_mat_oc_rc_filtered, 1, sd), decreasing = TRUE), 500)
expr_sub <- expr_mat_oc_rc_filtered[top_genes, ]

# Create annotation for columns
annotation_col <- data.frame(
  Group = sample_info$Group,
  Tissue = sample_info$Tissue,
  Status = sample_info$Status,
  SourceType = sample_info$SourceType
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
# check
any(is.na(sample_info$Group)) #FALSE
all(colnames(expr_mat_oc_rc_filtered) == rownames(sample_info)) #TRUE

# limma
design <- model.matrix(~0 + Group, data = sample_info)
colnames(design) <- levels(sample_info$Group)

colnames(design) <- make.names(colnames(design)) # for to make the name valid for limma contrasts TCGA-KIRC to TCGA.KIRC

contrast.matrix <- makeContrasts(OCCC_vs_ccRCC_ratio = (OCCC - GTEx_Ovary) - (TCGA.KIRC - GTEx_Renal_Cortex),
                                 levels = design)

fit <- lmFit(expr_mat_oc_rc_filtered, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend =TRUE)

DE_OCCC_vs_ccRCC_ratio <- topTable(fit2, coef="OCCC_vs_ccRCC_ratio", number=Inf, adjust.method="BH")
DE_OCCC_vs_ccRCC_ratio <- cbind(Geneid = rownames(DE_OCCC_vs_ccRCC_ratio), DE_OCCC_vs_ccRCC_ratio)
rownames(DE_OCCC_vs_ccRCC_ratio) <- NULL

# save DE results
wb <- createWorkbook()
addWorksheet(wb, "OCCC_vs_ccRCC_ratio")
writeData(wb, "OCCC_vs_ccRCC_ratio", DE_OCCC_vs_ccRCC_ratio)

write_csv(DE_OCCC_vs_ccRCC_ratio, file = "/Users/beyzaerkal/Desktop/occc_multi-omics/results/DE_results_OCCC_vs_ccRCC_ratio.csv")




####################################
# Only OCCC and GTEx Ovary samples
####################################


# using the already used expr_raw_oc_rc_filtered to subset only OCCC and GTEx Ovary samples for DE analysis -
# to avoid reprocessing the raw data again

# matrix
expr_mat_ovary <- expr_mat_oc_rc_filtered[, sample_info$Group %in% c("OCCC", "GTEx_Ovary")]

# metadata - and remove unused factors
sample_info_ovary <- droplevels(sample_info[sample_info$Group %in% c("OCCC", "GTEx_Ovary"), ])
all(colnames(expr_mat_ovary) == rownames(sample_info_ovary))

# save
#saveRDS(sample_info_ovary, file = "/Users/beyzaerkal/Desktop/occc_multi-omics/processed/sample_info_occc_ovary.rds")
#saveRDS(expr_mat_ovary, file = "/Users/beyzaerkal/Desktop/occc_multi-omics/processed/expression_matrix_occc_ovary.rds")

# model design : log2(OCCC / GTEx Ovary)
design_ovary <- model.matrix(~0 + Group, data = sample_info_ovary)
colnames(design_ovary) <- levels(sample_info_ovary$Group)
design_ovary

contrast.matrix_ovary <- makeContrasts(OCCC_vs_GTEx_Ovary_ratio = OCCC - GTEx_Ovary,
                                       levels = design_ovary)

fit <- lmFit(expr_mat_ovary, design_ovary)
fit2 <- contrasts.fit(fit, contrast.matrix_ovary)
fit2 <- eBayes(fit2, trend = TRUE)

# add the geneid to columns
DE_OCCC_vs_GTEx_Ovary_ratio <- topTable(fit2, coef="OCCC_vs_GTEx_Ovary_ratio", number=Inf, adjust.method="BH")
DE_OCCC_vs_GTEx_Ovary_ratio <- cbind(Geneid = rownames(DE_OCCC_vs_GTEx_Ovary_ratio), DE_OCCC_vs_GTEx_Ovary_ratio)
rownames(DE_OCCC_vs_GTEx_Ovary_ratio) <- NULL

# save
write_csv(DE_OCCC_vs_GTEx_Ovary_ratio, file = "/Users/beyzaerkal/Desktop/occc_multi-omics/results/DE_results_OCCC_vs_GTEx_Ovary_ratio.csv")


# save them as excel sheet
addWorksheet(wb, "OCCC_vs_GTEx_Ovary_ratio")
writeData(wb, "OCCC_vs_GTEx_Ovary_ratio", DE_OCCC_vs_GTEx_Ovary_ratio)

saveWorkbook(wb, "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/DE_results.xlsx", overwrite = TRUE)




##################
# immune cell infiltration save the OCCC expr mat in tsv format (linear TPM - without log2tr)
##################


#CCOC_df_mRNA
write_tsv(CCOC_df_mRNA, "/Users/beyzaerkal/Desktop/occc_multi-omics/inputs/CCOC_RNAseq_mRNA_TPM.tsv") # 200
hist(as.numeric(as.matrix(CCOC_df_mRNA[,-1])))
sum(duplicated(CCOC_df_mRNA$Geneid))

# ea_tpm_mRNA
hist(as.numeric(as.matrix(ea_tpm_mRNA[,-1])))
sum(duplicated(ea_tpm_mRNA$Geneid))
ea_tpm_mRNA <- ea_tpm_mRNA %>% group_by(Geneid) %>% summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
write_tsv(ea_tpm_mRNA, "/Users/beyzaerkal/Desktop/occc_multi-omics/inputs/E-MTAB_RNAseq_mRNA_TPM.tsv") # 7

# tpm_GSE160692_mRNA
write_tsv(tpm_GSE160692_mRNA, "/Users/beyzaerkal/Desktop/occc_multi-omics/inputs/GSE160692_RNA_seq_mRNA_TPM.tsv") # 11
hist(as.numeric(as.matrix(tpm_GSE160692_mRNA[,-1])))
sum(duplicated(tpm_GSE160692_mRNA$Geneid))

# tpm_GSE189553_mRNA
write_tsv(tpm_GSE189553_mRNA, "/Users/beyzaerkal/Desktop/occc_multi-omics/inputs/GSE189553_RNAseq_mRNA_TPM.tsv") # 11
hist(as.numeric(as.matrix(tpm_GSE189553_mRNA[,-1])))
sum(duplicated(tpm_GSE189553_mRNA$Geneid))

# tpm_SD1_OC_mRNA
write_tsv(tpm_SD1_OC_mRNA, "/Users/beyzaerkal/Desktop/occc_multi-omics/inputs/SD1_RNAseq_mRNA_TPM.tsv") # 6
hist(as.numeric(as.matrix(tpm_SD1_OC_mRNA[,-1])))
sum(duplicated(tpm_SD1_OC_mRNA$Geneid))












