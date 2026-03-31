
library(tidyverse)
library(readxl)
library(biomaRt)
library(purrr)
options(timeout = 300)
library(openxlsx)
library(tibble)
library(stringr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limma)
library(ggplot2)
library(impute)
library(preprocessCore)
set.seed(123)
###########################
# PROTEOMICS SECTION (OLD AND NEW PROTEOMICS DATA)
###########################


###########
# OLD
###########
# OCCC - proteomics

pro_L <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/raw_input_files/PXD032355_file/source_file_proL.xlsx", sheet = 14, col_names = TRUE)
# 40 samples in CCOC group, Protein expression, also have normal group (n=31)
#sheet 14 - raw protein intensity values from mass spectrometry PulseDIA
pro_L <- pro_L[, colSums(!is.na(pro_L)) > 0]
# metadata
meta <- dplyr::select(pro_L, X, pat_id, Histology8, Group, type)
unique(pro_L$type)
# filter 
pro_L_ccoc   <- pro_L %>% filter(type == "CCOC")
pro_L_normal <- pro_L %>% filter(type == "Normal")

# save 
#saveRDS(prot_ccoc, "proteomics_CCOC_only.rds")
#saveRDS(prot_normal, "proteomics_Normal_only.rds")
nrow(pro_L_ccoc)
nrow(pro_L_normal)

# remove metadata
pro_L_data <- dplyr::select(pro_L_ccoc, -pat_id, -Histology8, -Group, -type)

# pivot long and split uniprotid
pro_L_long <- pro_L_data %>% pivot_longer(-X, names_to = "Protein", values_to = "Intensity") %>%
  separate(Protein, into = c("UniProtID", "Gene"), sep = "_")

# get gene level by mean - avoid duplicates
pro_L_gene <- pro_L_long %>% group_by(Gene, X) %>%
  summarise(Intensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

# pivot back to wide
pro_L_matrix <- pro_L_gene %>%
  pivot_wider(names_from = X, values_from = Intensity) %>%
  column_to_rownames("Gene")
# gene names as a column
pro_L_matrix_out <- pro_L_matrix %>%
  as.data.frame(check.names = FALSE) %>%
  tibble::rownames_to_column("Geneid")

wb <- createWorkbook()
addWorksheet(wb, "Qian et al. (2024)")
writeData(wb, "Qian et al. (2024)", pro_L_matrix_out)
# 5 genes only, 40 samples
#write.csv(pro_L_matrix_out, "/Users/beyzaerkal/Desktop/internship/internship_env/Input_OC_prot/Proteomic_land_PROT.csv", row.names = FALSE)


##################################
# PROT: DepMap OCCC cell lines
##################################

# only have RPPA 
prot_ALL <- read_csv("/Users/beyzaerkal/Desktop/internship/internship_env/raw_input_files/DepMap_Prot/Harmonized_RPPA_CCLE_subsetted.csv", col_names = TRUE)

OCids <- c("ACH-000885", "ACH-000906", "ACH-000719", "ACH-000527", "ACH-000663", "ACH-000324", "ACH-000646")

prot_OCCC <- prot_ALL[prot_ALL$depmap_id %in% OCids, ]

# remove extra col
prot_OCCC_clean <- prot_OCCC[ , !colnames(prot_OCCC) %in% c("depmap_id", grep("lineage", colnames(prot_OCCC), value = TRUE))]
# pivot and prot id+delete colname
long <- prot_OCCC_clean %>% pivot_longer(cols = -cell_line_display_name, names_to = "ProtID", values_to = "value")

# cell line -> colnames
prot_matrix <- long %>% pivot_wider(names_from = cell_line_display_name, values_from = value)

prot_matrix <- as.data.frame(prot_matrix)
rownames(prot_matrix) <- prot_matrix$ProtID
prot_matrix$ProtID <- NULL


# get UniprotID to gene symbol
prot_matrix$uniprot_id <- sub(" .*", "", rownames(prot_matrix))

length(prot_matrix$uniprot_id) # 144
# get gene sym
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

annot <- getBM(attributes = c("uniprotswissprot", "hgnc_symbol", "gene_biotype"),
               filters = "uniprotswissprot",
               values = prot_matrix$uniprot_id,
               mart = ensembl)

annot_pc <- annot[annot$gene_biotype == "protein_coding", ]

#merge into my data
merged <- merge(prot_matrix, annot_pc[, c("uniprotswissprot", "hgnc_symbol")],by.x = "uniprot_id", by.y = "uniprotswissprot", all.x = FALSE)

rownames(merged) <- merged$hgnc_symbol
# drop
merged$hgnc_symbol <- NULL
merged$uniprot_id <- NULL

write.csv(merged, "DepMap_hRPPA_CCLE_OCCC_filtered.csv") # 144

addWorksheet(wb, "Nusinow et al. (2020)")
writeData(wb, "Nusinow et al. (2020)", merged)


##################
# ccRCC proteomics
##################

prot_renal_raw <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/raw_files_ccRCC/TCGA-KIRC.protein.tsv")

prot_renal_raw[is.na(prot_renal_raw)] <- 0

prot_renal_raw$CleanName <- str_replace(prot_renal_raw$peptide_target, "[-_].*$", "")

# mapping protein symbols using ALIAS2EG and then back to SYMBOL
mapped <- mapIds(org.Hs.eg.db, 
                 keys = prot_renal_raw$CleanName,
                 column = "SYMBOL",
                 keytype = "ALIAS",
                 multiVals = "first")

prot_renal_raw$GeneSymbol <- mapped[prot_renal_raw$CleanName] # 487

unmapped <- prot_renal_raw[is.na(prot_renal_raw$GeneSymbol), "peptide_target"]
length(unmapped)
head(unmapped) #161

prot_renal_mapped <- prot_renal_raw[!is.na(prot_renal_raw$GeneSymbol), ] # 326


colnames(prot_renal_mapped)[colnames(prot_renal_mapped) == "peptide_target"] <- "Geneid"
prot_renal_mapped$Geneid <- toupper(prot_renal_mapped$Geneid) # for capital cases

# clean the Geneid column 
prot_renal_mapped$Geneid <- toupper(gsub("[-_].*$", "", prot_renal_mapped$Geneid))

write.csv(prot_renal_mapped, "TCGA-KIRC.protein_with_genes2.csv", row.names = FALSE)

addWorksheet(wb, "TCGA-KIRC_RPPA")
writeData(wb, "TCGA-KIRC_RPPA", prot_renal_mapped)






##################
#NEW AND USED DATASETS BELOW
##################





###################
# The proteome of CCOC paper PRIDE ARCHIVE Protein txt files with different batches
###################


path <- "/Users/beyzaerkal/Desktop/internship/internship_env/raw_input_files/prot_OCCC_JJ2022_raw"

files <- list.files(path, full.names = TRUE)

read_oc_p <- function(file) {
  df <- read.delim(file, sep = "\t", quote = "\"", check.names = FALSE)
  
  colnames(df) <- gsub('^\"|\"$', '', colnames(df)) # clean col
  abundance_cols <- grep("^Abundance", names(df), value = TRUE) # abundance col
  
  df_long <- df %>% dplyr::select(Accession, all_of(abundance_cols)) %>%
    pivot_longer(cols = all_of(abundance_cols), names_to = "Channel", values_to = "Abundance") %>%
    mutate(Abundance = trimws(Abundance),
           Abundance = na_if(Abundance, ""), 
           Abundance = na_if(Abundance, "Not Found"),
           Abundance = as.numeric(Abundance), 
           #Abundance = ifelse(is.na(Abundance), 0, Abundance),
           Batch = stringr::str_extract(Channel, "F[0-9]+"),
           Label = stringr::str_extract(Channel, "12[6-9]|13[0-1][NC]?"),
           sample_id = paste0(Batch, "_", Label))
  
  return(df_long)
}


# all in one df 
proteomics_long <- purrr::map_df(files, read_oc_p)
any(is.na(proteomics_long))
# standardise - remove N/C
proteomics_long <- proteomics_long %>% mutate(sample_id = gsub("(N|C)$", "", sample_id))

# plot density plot for checks
ggplot(proteomics_long, aes(x = log2(Abundance + 1), group = sample_id, colour = sample_id)) +
  geom_density(alpha = 0.3) + theme_bw() +
  labs(x = "log2 Abundance", y = "Density") + theme(legend.position = "none")

# aggregate multiple channels for same sample
proteomics_summary <- proteomics_long %>%
  group_by(sample_id, Accession) %>%
  summarise(Abundance = median(Abundance, na.rm = TRUE), .groups = "drop")

# pivot wide
prot_matrix <- proteomics_summary %>%
  pivot_wider(names_from = Accession, values_from = Abundance, values_fill = NA) 

# log2 transform
#prot_log2 <- prot_matrix %>% mutate(across(-sample_id, ~ log2(.x + 1)))
prot_log2 <- prot_matrix %>%mutate(across(-sample_id, ~ log2(as.numeric(.x)))) # keep NA

######## check TMT batch effects with PCA
pca_matrix <- prot_log2 %>%
  tibble::column_to_rownames("sample_id")

# numeric
pca_matrix <- data.matrix(pca_matrix)

# remove zero variance
pca_matrix <- pca_matrix[, apply(pca_matrix, 2, sd, na.rm = TRUE) > 0]

# remove NA or Inf
pca_matrix <- pca_matrix[, colSums(is.na(pca_matrix) | is.infinite(pca_matrix)) == 0]

pca <- prcomp(pca_matrix, scale. = TRUE)

pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                     sample_id = rownames(pca$x))

# extract batch
pca_df$Batch <- stringr::str_extract(pca_df$sample_id, "F[0-9]+")
# plot PCA
ggplot(pca_df, aes(PC1, PC2, color = Batch)) + 
  geom_point(size = 3) +
  theme_bw() +
  labs(x = "PC1", y = "PC2")

###############

# uniprot to gene symbol
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
listAttributes(mart)[grep("uniprot", listAttributes(mart)$name), ]

uniprot_ids <- unique(proteomics_long$Accession)
# Swiss-prot
mapping_swiss <- getBM(attributes = c("uniprotswissprot", "hgnc_symbol"),
                       filters = "uniprotswissprot", values = uniprot_ids, mart = mart)
#head(mapping_swiss)
#TrEMBL
mapping_trembl <- getBM(attributes = c("uniprotsptrembl"), values = uniprot_ids, mart = mart)
#head(mapping_trembl)

mapping_trembl <- mapping_trembl %>% dplyr::rename(Accession = uniprotsptrembl) %>% mutate(SYMBOL = NA_character_)
# combine
mapping <- dplyr::bind_rows(mapping_swiss %>% dplyr::rename(Accession = uniprotswissprot, SYMBOL = hgnc_symbol),
                            mapping_trembl)

# keep rows with SYMBOL first (Swiss-Prot)
mapping <- mapping %>% arrange(!is.na(SYMBOL)) %>% distinct(Accession, .keep_all = TRUE)

# convert accession matrix to gene matrix 
prot_long <- prot_log2 %>%
  pivot_longer(cols = -sample_id, names_to = "Accession", values_to = "Abundance") %>%
  group_by(sample_id, Accession) %>%
  summarise(Abundance = median(Abundance, na.rm = TRUE), .groups = "drop") %>% 
  inner_join(mapping %>% filter(!is.na(SYMBOL)), by = "Accession")

# aggregate multi channel samples per protein using median
prot_gene <- prot_long %>%
  group_by(sample_id, SYMBOL) %>%
  summarise(Abundance = median(Abundance, na.rm = TRUE), .groups = "drop") %>% 
  filter(!is.na(SYMBOL) & SYMBOL != "") %>%
  pivot_wider(names_from = sample_id, values_from = Abundance) %>% # , values_fill = 0
  column_to_rownames(var = "SYMBOL") 

prot_gene <- prot_gene %>% tibble::rownames_to_column(var = "Geneid")

# convert to matrix for centering
prot_mat <- prot_gene %>% tibble::column_to_rownames("Geneid") %>% as.matrix()

# filter low variance
keep <- rowMeans(!is.na(prot_mat)) >= 0.7
prot_mat <- prot_mat[keep, ]

# impute missing values
prot_mat <- impute.knn(prot_mat)$data

# adding median centering for consistent data types
# median-center per protein/gene not sample
prot_centered <- sweep(prot_mat, 1, apply(prot_mat, 1, median, na.rm = TRUE), "-")

# convert back to df
prot_gene_centered <- as.data.frame(prot_centered) %>%
  tibble::rownames_to_column("Geneid")

# save both log2 and centered
write_csv(prot_gene, "/Users/beyzaerkal/Desktop/occc_multi-omics/processed/proteomics_JJ2022_log2.csv")
write_csv(prot_gene_centered, "/Users/beyzaerkal/Desktop/occc_multi-omics/processed/proteomics_JJ2022_log2_centered.csv")

# svae to excel
wb <- createWorkbook()
addWorksheet(wb, "Ji et al. (2022) log2")
writeData(wb, "Ji et al. (2022) log2", prot_gene)

addWorksheet(wb, "Ji et al. (2022) centered")
writeData(wb, "Ji et al. (2022) centered", prot_gene_centered)


##################
# ccRCC prot
# PDC000127  - tmt10 - protoemics data commons
####################

ccrcc_prot_raw <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/raw_files_ccRCC/CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv")

ccrcc_prot <- ccrcc_prot_raw %>% filter(!Gene %in% c("Mean", "Median", "StdDev")) # remove summary rows

# only log ratio cols
ccrcc_prot <- ccrcc_prot %>% dplyr::select(Gene, contains("Log Ratio")) %>% dplyr::select(-contains("Unshared"))

ccrcc_matrix <- ccrcc_prot %>% column_to_rownames("Gene") %>% as.matrix()
# filter
keep <- rowMeans(!is.na(ccrcc_matrix)) >= 0.7
ccrcc_matrix <- ccrcc_matrix[keep, ]
# impute missing values
ccrcc_matrix <- impute.knn(ccrcc_matrix)$data

any(is.na(ccrcc_matrix))

# check distribution
boxplot(ccrcc_matrix, las = 2, main = "ccRCC CPTAC log ratios")

# median-center ccRCC samples 
ccrcc_centered <- sweep(ccrcc_matrix, 1, apply(ccrcc_matrix, 1, median, na.rm = TRUE), "-")

ccrcc_centered_df <- as.data.frame(ccrcc_centered) %>% rownames_to_column("Geneid")

# already HGNC gene symbols
# save both log2 ratio and centered have smae values meaning it was already cnetered around median (ccrcc_prot vs ccrcc_prot_centered)
#write_csv(ccrcc_prot, "/Users/beyzaerkal/Desktop/occc_multi-omics/processed/CPTAC_ccRCC_proteomics_log_ratios_raw.csv")
write_csv(ccrcc_centered_df, "/Users/beyzaerkal/Desktop/occc_multi-omics/processed/CPTAC_ccRCC_proteomics_log_ratios_centered.csv")

#addWorksheet(wb, "CPTAC ccRCC")
#writeData(wb, "CPTAC ccRCC", ccrcc_prot)

addWorksheet(wb, "CPTAC ccRCC centered")
writeData(wb, "CPTAC ccRCC centered", ccrcc_centered_df)

#saveWorkbook(wb, file = "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary", overwrite = TRUE)

######################
# renormalise for merging - otherwise causes negative domination - check on hist
# DE limma
######################
occc_mat <- prot_gene_centered %>% column_to_rownames("Geneid") %>% as.matrix()

ccrcc_mat <- ccrcc_centered_df %>% column_to_rownames("Geneid") %>% as.matrix()

common_genes <- intersect(rownames(occc_mat), rownames(ccrcc_mat))

occc_mat <- occc_mat[common_genes, ]
ccrcc_mat <- ccrcc_mat[common_genes, ]

prot_combined <- cbind(occc_mat, ccrcc_mat)
prot_combined_qn <- normalize.quantiles(as.matrix(prot_combined))

rownames(prot_combined_qn) <- rownames(prot_combined)
colnames(prot_combined_qn) <- colnames(prot_combined)

group <- factor(c(rep("OCCC", ncol(occc_mat)), rep("ccRCC", ncol(ccrcc_mat)))) 
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(prot_combined_qn, design)
contrast.matrix <- makeContrasts(OCCC_vs_ccRCC = OCCC - ccRCC, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, coef = "OCCC_vs_ccRCC", number = Inf, adjust.method = "BH")
results_df <- results %>% tibble::rownames_to_column(var = "Geneid")

write_csv(results_df, "/Users/beyzaerkal/Desktop/occc_multi-omics/results/proteomics_results/DE_OCCC_vs_ccRCC_proteomics_qn.csv")
hist(results$logFC, breaks = 50)

addWorksheet(wb, "Combined_matrix")
writeData(wb, "Combined_matrix", prot_combined_qn_df)
saveWorkbook(wb, file = "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary", overwrite = TRUE)


prot_combined_qn_df <- as.data.frame(prot_combined_qn) %>%
  tibble::rownames_to_column(var = "Geneid")

# FOR DE
wb2 <- createWorkbook()
addWorksheet(wb2, "DE_prot_OCCC_vs_ccRCC")
writeData(wb2, "DE_prot_OCCC_vs_ccRCC", results_df)
saveWorkbook(wb2, file = "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary", overwrite = TRUE)


# checking the scaling for both ccRCC and OCCC
par(mfrow=c(1,2))
# OCCC raw log2 (no median centering)
boxplot(as.matrix(prot_gene %>% column_to_rownames("Geneid")), main="OCCC raw log2", las=2, cex.axis=0.7)
# ccRCC raw log ratios - no median cnetreing 
boxplot(ccrcc_matrix, main="ccRCC log ratios", las=2, cex.axis=0.7)
#############
# median centreboard
par(mfrow=c(1,2))
boxplot(as.matrix(prot_gene_centered %>% column_to_rownames("Geneid")), main="OCCC log2 median centered", las=2, cex.axis=0.7)
boxplot(ccrcc_centered, main="ccRCC log ratios median centered", las=2, cex.axis=0.7)




