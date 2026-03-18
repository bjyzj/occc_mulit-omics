

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(limma)
library(readxl)
library(biomaRt)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(DOSE)
set.seed(123)

## Load data
occc_prot <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/proteomics_JJ2022_log2_centered.csv")
ccrcc_prot <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/CPTAC_ccRCC_proteomics_log_ratios_centered.csv")

# send genes to rownames
occc_mat <- occc_prot %>% column_to_rownames("Geneid") %>% as.matrix()
ccrcc_mat <- ccrcc_prot %>% column_to_rownames("Geneid") %>% as.matrix()

# expression matrix and intersection of genes
common_genes <- intersect(rownames(occc_mat), rownames(ccrcc_mat))
occc_mat <- occc_mat[common_genes, ]
ccrcc_mat <- ccrcc_mat[common_genes, ]

prot_combined <- cbind(occc_mat, ccrcc_mat)

#label
group <- c(rep("OCCC", ncol(occc_mat)), rep("ccRCC", ncol(ccrcc_mat)))
group <- factor(group)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

table(group)
#limma
fit <- lmFit(prot_combined, design)

# contrast matrix
contrast.matrix <- makeContrasts(OCCC_vs_ccRCC = OCCC - ccRCC, levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, coef = "OCCC_vs_ccRCC", number = Inf, adjust.method = "BH")

#write.csv(results, "/Users/beyzaerkal/Desktop/occc_multi-omics/results/proteomics_results/DE_OCCC_vs_ccRCC_proteomics.csv")


###########################
# ORA
###########################
results <- results %>% tibble::rownames_to_column(var = "Gene")
gene_map <- bitr(results$Gene,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

results_mapped <- results %>% inner_join(gene_map, by = c("Gene" = "SYMBOL"))

up_genes <- results_mapped %>% filter(adj.P.Val < 0.05 & logFC > 1)
down_genes <- results_mapped %>% filter(adj.P.Val < 0.05 & logFC < -1)

up_entrez_ids <- unique(up_genes$ENTREZID)
down_entrez_ids <- unique(down_genes$ENTREZID)


ekegg_up <- enrichKEGG(gene = up_entrez_ids,
                       organism = "hsa",
                       pvalueCutoff = 0.05)

reactome_up <- enrichPathway(gene = up_entrez_ids,
                             organism = "human",
                             pvalueCutoff = 0.05)



ekegg_down <- enrichKEGG(gene = down_entrez_ids,
                         organism = "hsa",
                         pvalueCutoff = 0.05)

reactome_down <- enrichPathway(gene = down_entrez_ids,
                               organism = "human",
                               pvalueCutoff = 0.05)

dotplot(ekegg_up, showCategory = 20) + theme_minimal()
dotplot(reactome_up, showCategory = 20) + theme_minimal()

dotplot(ekegg_down, showCategory = 20) + theme_minimal()
dotplot(reactome_down, showCategory = 20) + theme_minimal()


################
# GSEA
################

# ranked gene list for GSEA
all_genes <- results_mapped
gene_list <- all_genes$logFC
names(gene_list) <- all_genes$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)  

gene_list <- results_mapped$t  # t stats
names(gene_list) <- results_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

sum(duplicated(names(gene_list))) 
gene_list <- gene_list[!duplicated(names(gene_list))]  
# GO BP GSEA
gsea_go <- gseGO(geneList = gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

# KEGG GSEA
gsea_kegg <- gseKEGG(geneList = gene_list,
                     organism = "hsa",
                     minGSSize = 10,
                     pvalueCutoff = 0.05)

# Reactome GSEA
gsea_reactome <- gsePathway(geneList = gene_list,
                            organism = "human",
                            pvalueCutoff = 0.05,
                            verbose = FALSE)

dotplot(gsea_go, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_kegg, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_reactome, showCategory = 20, color = "NES") + theme_minimal()

#require(DOSE)
dotplot(gsea_reactome, showCategory=10, split=".sign") + facet_grid(.~.sign)


#####################
#kinases
####################

# kinase list ( threshold 0.05)
kinases <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/kinase_basic.xlsx", col_names = TRUE)
colnames(kinases)

# SYMBOL col match kinase list
up_genes$SYMBOL <- toupper(up_genes$Gene)
down_genes$SYMBOL <- toupper(down_genes$Gene)

kinases$Offical_gene_symbol <- toupper(kinases$Offical_gene_symbol)

# filter for kinases
up_kinases <- up_genes %>% dplyr::filter(SYMBOL %in% kinases$Offical_gene_symbol)
down_kinases <- down_genes %>% dplyr::filter(SYMBOL %in% kinases$Offical_gene_symbol)

# ENTREZ IDs
up_kinases_ids <- unique(up_kinases$ENTREZID)
down_kinases_ids <- unique(down_kinases$ENTREZID)


kegg_up_kinase <- enrichKEGG(gene = up_kinases_ids,
                             organism = "hsa",
                             pvalueCutoff = 0.05)

kegg_down_kinase <- enrichKEGG(gene = down_kinases_ids,
                               organism = "hsa",
                               pvalueCutoff = 0.05)

reactome_up_kinase <- enrichPathway(gene = up_kinases_ids,
                                    organism = "human",
                                    pvalueCutoff = 0.05)

reactome_down_kinase <- enrichPathway(gene = down_kinases_ids,
                                      organism = "human",
                                      pvalueCutoff = 0.05)


dotplot(kegg_up_kinase, showCategory = 20) + theme_minimal()
dotplot(reactome_up_kinase, showCategory = 20) + theme_minimal()

dotplot(kegg_down_kinase, showCategory = 20) + theme_minimal()
dotplot(reactome_down_kinase, showCategory = 20) + theme_minimal()

