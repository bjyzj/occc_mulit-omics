
library(tidyverse)
library(readxl)
library(biomaRt)
library(openxlsx)
library(SummarizedExperiment)
library(GenomicRanges)
library(tibble)
library(stringr)
library(org.Hs.eg.db)
library(AnnotationDbi)

#####################################
# Molecular Subclasses of Clear Cell Ovarian Carcinoma and Their Impact on Disease Behavior and Outcomes
# (1) sample 211 primary - tpm
#####################################

CCOC_rds_file <- readRDS("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/CCOC.RNAseq.normalized.rds")

# remove _PAR_Y
head(rownames(CCOC_rds_file))
rownames(CCOC_rds_file) <- sub("_PAR_Y$", "", rownames(CCOC_rds_file))
# remove verison
rds_version_ids <- rownames(CCOC_rds_file)
rds_gene_ids <- sub("\\..*$", "", rds_version_ids)
# all are ensg -> gene symbol (onyl mrna)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

bm_results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
                    filters = "ensembl_gene_id", values = rds_gene_ids, mart = ensembl)

bm_results_mrna <- bm_results %>% filter(gene_biotype == "protein_coding")
bm_results_mrna <- bm_results_mrna %>% filter(hgnc_symbol != "" & !is.na(hgnc_symbol))

# merge back
mapping_df <- data.frame(ensembl_with_version = rds_version_ids, ensembl_gene_ids = rds_gene_ids)
# only protein-coding genes to mapping_df
colnames(mapping_df)[colnames(mapping_df) == "ensembl_gene_ids"] <- "ensembl_gene_id"
mapping_df <- merge(mapping_df, bm_results_mrna, by = "ensembl_gene_id",all.x = FALSE)

CCOC_rds_file <- CCOC_rds_file[mapping_df$ensembl_with_version, ]

CCOC_df <- as.data.frame(CCOC_rds_file)
# duplicates

CCOC_df$hgnc_symbol <- mapping_df$hgnc_symbol
CCOC_df <- CCOC_df[, c(ncol(CCOC_df), 1:(ncol(CCOC_df)-1))]
# take mean of duplicates- mean for each numeric column separately
CCOC_df <- CCOC_df %>% group_by(hgnc_symbol) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% ungroup()

CCOC_df <- as.data.frame(CCOC_df)
# nrow = 19,330
# chaneg col name
colnames(CCOC_df) <- trimws(colnames(CCOC_df))
colnames(CCOC_df)[colnames(CCOC_df) == "hgnc_symbol"] <- "Geneid"

write.csv(CCOC_df,
          file = "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/CCOC_RNAseq_mRNA_TPM.csv",
          row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "Bolton et al., 2022")
writeData(wb, "Bolton et al., 2022", CCOC_df)

##########################################
# (2) E-MTAB (expression atlas - occc and their types) TPM
##########################################

ea_tpm <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/E-MTAB-2770-query-results.tsv", skip = 4, col_names = TRUE)

# NAs fill with 0 - assume no expression 
ea_tpm[ , 3:ncol(ea_tpm)][is.na(ea_tpm[ , 3:ncol(ea_tpm)])] <- 0
# skip first 2 columns as its id only
ea_tpm <- ea_tpm[, -1]
colnames(ea_tpm) <- gsub(" ", "_", colnames(ea_tpm))
colnames(ea_tpm)[colnames(ea_tpm) == "Gene_Name"] <- "Geneid"

sum(duplicated(ea_tpm)) #132
ea_tpm_unique <- ea_tpm %>% group_by(Geneid) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% ungroup()
# mRNA filter
genes <- unique(ea_tpm$Geneid)
annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = genes,
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
ea_tpm_mRNA <- ea_tpm %>% filter(Geneid %in% mRNA_genes)


write.csv(ea_tpm_mRNA,
          file = "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/E-MTAB_RNAseq_mRNA_TPM.csv",
          row.names = FALSE)

addWorksheet(wb, "E-MTAB-2770")
writeData(wb, "E-MTAB-2770", ea_tpm_mRNA)

#######################################
# (3) GSE160692_data- raw to tpm conv
#######################################

GSE160692_data <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE160692_all/GSE160692_OVA_UTE_Raw_Transcripts_GEO.tsv")
table(GSE160692_data$`Gene Symbol`)
duplicate_gene_GSE160692 <- GSE160692_data$`Gene Symbol`[duplicated(GSE160692_data$`Gene Symbol`)]

# only the OVA columns
ova_cols <- grep("^OVA", names(GSE160692_data), value = TRUE)

# group by all non-OVA columns, and average only OVA columns
duplicate_means_GSE160692 <- GSE160692_data %>%
  group_by(across(-all_of(ova_cols))) %>%
  summarise(across(all_of(ova_cols), mean), .groups = "drop")

ova_table <- dplyr::select(duplicate_means_GSE160692, `Gene Symbol`, Start, End, Chromosome, Strand, dplyr::all_of(ova_cols))

# there are gene symbols with dates that start with number might not match with others
# during alignment only align the ones match together

gene_with_dates <- ova_table %>% filter(grepl("^[0-9]", `Gene Symbol`)) # viewed first

GSE160692_raw_count <- ova_table %>% filter(!grepl("^[0-9]", `Gene Symbol`))

# raw to tpm
gene_id = GSE160692_raw_count$`Gene Symbol`

gr <- GRanges(seqnames = GSE160692_raw_count$Chromosome,
              ranges = IRanges(start = GSE160692_raw_count$Start, end = GSE160692_raw_count$End),
              strand = GSE160692_raw_count$Strand,
              gene_id = GSE160692_raw_count$`Gene Symbol`)


se <- SummarizedExperiment(assays = list(counts = as.matrix(GSE160692_raw_count[, ova_cols])),
                           rowRanges = gr,
                           colData = DataFrame(samples = ova_cols))

counts <- assay(se, "counts")
lengths_kb <- width(rowRanges(se)) / 1000

count2tpm <- function(counts, lengths_kb) {
  rpk <- counts / lengths_kb
  tpm <- t( t(rpk) / colSums(rpk) ) * 1e6
  return(tpm)
}

tpm_GSE160692 <- count2tpm(counts, lengths_kb)
rownames(tpm_GSE160692) <- GSE160692_raw_count$`Gene Symbol`

tpm_GSE160692 <- as.data.frame(tpm_GSE160692)
tpm_GSE160692<- cbind(`Gene Symbol` = rownames(tpm_GSE160692), tpm_GSE160692)
rownames(tpm_GSE160692) <- NULL

colnames(tpm_GSE160692) <- trimws(colnames(tpm_GSE160692))
colnames(tpm_GSE160692)[grep("Gene", colnames(tpm_GSE160692))] <- "Geneid"
# mRNA filter
annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = unique(tpm_GSE160692$Geneid),
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
tpm_GSE160692_mRNA <- tpm_GSE160692 %>% filter(Geneid %in% mRNA_genes)

write.csv(tpm_GSE160692_mRNA, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE160692_RNA_seq_mRNA_TPM.csv", row.names = FALSE)

addWorksheet(wb, "GSE160692")
writeData(wb, "GSE160692", tpm_GSE160692_mRNA)

########################################
# (4) tpm file extract ccc GSE189553
########################################

GSE189553_data <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE189553_all/GSE189553_gene_TPM_matrix.txt")

GSE189553_data <- dplyr::select(GSE189553_data, starts_with("gene"), starts_with("CCC"))

GSE189553_data <- dplyr::rename(GSE189553_data, Geneid = gene) 
GSE189553_data<- dplyr::relocate(GSE189553_data, Geneid, .before = 1)

tpm_GSE189553 <- as.data.frame(GSE189553_data)

sum(duplicated(tpm_GSE189553$Geneid)) # 1706

tpm_GSE189553 <- tpm_GSE189553 %>% group_by(Geneid) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

sum(duplicated(tpm_GSE189553$Geneid))

tpm_GSE189553 <- tpm_GSE189553 %>% filter(!grepl("^[0-9]", Geneid))
# mRNA filter
annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = unique(tpm_GSE189553$Geneid),
               mart = ensembl)


mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
tpm_GSE189553_mRNA <- tpm_GSE189553 %>% filter(Geneid %in% mRNA_genes)

write.csv(tpm_GSE189553_mRNA, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/GSE189553_RNAseq_mRNA_TPM.csv", row.names = FALSE)

addWorksheet(wb, "GSE189553")
writeData(wb, "GSE189553", tpm_GSE189553_mRNA)

###################################
# (5) rpkm to tpm (113 sample)
###################################

rpkm_113_data <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/113_OCCC/RPKM_expdata_113OCCC.xlsx", sheet = 1, col_names = TRUE)

colnames(rpkm_113_data) <- trimws(colnames(rpkm_113_data)) # colname had issue so cleaned
colnames(rpkm_113_data)[1] <- "Geneid"

rpkm_matrix <- as.matrix(rpkm_113_data[, -1])

rpkm2tpm <- function(rpkm_matrix) {
  tpm <- apply(rpkm_matrix, 2, function(x) {
    x / sum(x, na.rm = TRUE) * 1e6})
  return(tpm)}

tpm_matrix <- rpkm2tpm(rpkm_matrix)
tpm_113 <- cbind(Geneid = rpkm_113_data$Geneid, as.data.frame(tpm_matrix))
# no NAs 

sum(duplicated(tpm_113$Geneid))
# mRNA filter
annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = unique(tpm_113$Geneid),
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)

tpm_113_mRNA <- tpm_113 %>% filter(Geneid %in% mRNA_genes) 

write.csv(tpm_113_mRNA, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/tpm113_RNAseq_mRNA_TPM.csv", row.names = FALSE)

addWorksheet(wb, "Stružinská et al., 2023")
writeData(wb, "Stružinská et al., 2023", tpm_113_mRNA)

####################################
# (6) rpkm to tpm # paper: Systematic Identification of Characteristic Genes of Ovarian Clear Cell Carcinoma
####################################

rpkm_SD1 <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/ijms-20-04330-s001/Supplmentary Data S1.xlsx", sheet = 1, skip = 2)

rpkm_SD1_OC <- dplyr::select(rpkm_SD1, starts_with("symbol"), starts_with("C"))
rpkm_SD1_Norm <- dplyr::select(rpkm_SD1, starts_with("symbol"), starts_with("N"))

Geneid <- rpkm_SD1_OC[[1]]
rpkm_matrix <- as.matrix(rpkm_SD1_OC[, -1])


rpkm2tpm <- function(rpkm) {
  tpm <- apply(rpkm_matrix, 2, function(x) x / sum(x, na.rm = TRUE) * 1e6)
  return(tpm)}

tpm_matrix <- rpkm2tpm(rpkm_matrix)
tpm_SD1_OC <- cbind(Geneid = Geneid, as.data.frame(tpm_matrix))
tpm_SD1_OC <- as.data.frame(tpm_SD1_OC)

#duplicates
tpm_SD1_OC <- tpm_SD1_OC %>% group_by(Geneid) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

sum(duplicated(tpm_SD1_OC$Geneid))
# mRNA filter
annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = unique(tpm_SD1_OC$Geneid),
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
tpm_SD1_OC_mRNA <- tpm_SD1_OC %>% filter(Geneid %in% mRNA_genes) 

write.csv(tpm_SD1_OC_mRNA, "/Users/beyzaerkal/Desktop/internship/internship_env/OCCC_datasets/SD1_RNAseq_mRNA_TPM.csv", row.names = FALSE)

addWorksheet(wb, "Nagasawa et al., 2019")
writeData(wb, "Nagasawa et al., 2019", tpm_SD1_OC_mRNA)

###################################
# # Kidney
###################################
tpm_ccRCC <- read_tsv("/Users/beyzaerkal/Desktop/internship/internship_env/raw_files_ccRCC/TCGA-KIRC.star_tpm.tsv")

tpm_ccRCC$Ensembl_ID <- sub("\\..*", "", tpm_ccRCC$Ensembl_ID)

sum(duplicated(tpm_ccRCC$Ensembl_ID))

tpm_ccRCC_mean <- tpm_ccRCC %>% group_by(Ensembl_ID) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_map <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "ensembl_gene_id",
                   values = tpm_ccRCC_mean$Ensembl_ID,
                   mart = ensembl)

tpm_ccRCC <- tpm_ccRCC_mean %>% left_join(gene_map, by = c("Ensembl_ID" = "ensembl_gene_id"))

# gene symbols have NAs
sum(is.na(tpm_ccRCC)) # [1] 1253
sum(is.na(tpm_ccRCC$hgnc_symbol)) # [1] 1253

# changes to columns
tpm_ccRCC$Geneid <- tpm_ccRCC$hgnc_symbol
tpm_ccRCC$hgnc_symbol <- NULL
tpm_ccRCC$Ensembl_ID <- NULL
tpm_ccRCC <- tpm_ccRCC %>% dplyr:: select(Geneid, everything())

tpm_ccRCC <- as.data.frame(tpm_ccRCC)

#duplicates
tpm_ccRCC <- tpm_ccRCC %>% group_by(Geneid) %>% 
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

sum(duplicated(tpm_ccRCC$Geneid))

annot <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "hgnc_symbol",
               values = unique(tpm_ccRCC$Geneid),
               mart = ensembl)

mRNA_genes <- annot %>% filter(gene_biotype == "protein_coding") %>% pull(hgnc_symbol)
# filter 
tpm_ccRCC_mRNA <- tpm_ccRCC %>% filter(Geneid %in% mRNA_genes) 

write.csv(tpm_ccRCC_mRNA, "/Users/beyzaerkal/Desktop/internship/internship_env/input_ccRCC/ccRCC_RNAseq_mRNA_TPM.csv", row.names = FALSE)

addWorksheet(wb, "TCGA-KIRC")
writeData(wb, "TCGA-KIRC", tpm_ccRCC_mRNA)

saveWorkbook(wb, file = "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary", overwrite = TRUE)






