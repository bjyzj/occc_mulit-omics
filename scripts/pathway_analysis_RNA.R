
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(readxl)
set.seed(123)
# Load the data
DE_OCvsRC <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/transcriptomics_results/DE_results_OCCC_vs_ccRCC_ratio.csv")
DE_OCvsGTExOV <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/transcriptomics_results/DE_results_OCCC_vs_GTEx_Ovary_ratio.csv")

#####################
# OC vs RC ORA

DE_OCvsRC$SYMBOL <- DE_OCvsRC$Geneid

up_genes <- DE_OCvsRC[DE_OCvsRC$adj.P.Val < 0.05 & DE_OCvsRC$logFC > 1,]
down_genes <- DE_OCvsRC[DE_OCvsRC$adj.P.Val < 0.05 & DE_OCvsRC$logFC < -1,]

up_entrez <- bitr(up_genes$SYMBOL,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)

down_entrez <- bitr(down_genes$SYMBOL,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

# remove duplicates
sum(duplicated(up_genes$SYMBOL)) #0
sum(duplicated(down_genes$SYMBOL))#0
# remove duplicates in case
up_entrez_ids <- unique(up_entrez$ENTREZID)
down_entrez_ids <- unique(down_entrez$ENTREZID)

###
kegg_up <- enrichKEGG(gene = up_entrez_ids,
                      organism = "hsa",
                      pvalueCutoff = 0.05)

kegg_down <- enrichKEGG(gene = down_entrez_ids,
                        organism = "hsa",
                        pvalueCutoff = 0.05)
###
reactome_up <- enrichPathway(gene = up_entrez_ids,
                             organism = "human",
                             pvalueCutoff = 0.05)

reactome_down <- enrichPathway(gene = down_entrez_ids,
                               organism = "human",
                               pvalueCutoff = 0.05)
##
dotplot(kegg_up, showCategory = 20)  + theme_minimal()
dotplot(reactome_up, showCategory = 20) + theme_minimal()
                
dotplot(kegg_down, showCategory = 20) + theme_minimal()
dotplot(reactome_down, showCategory = 20)  + theme_minimal()

                                
#####################
# OC vs RC GSEA
# ranked gene list for GSEA

all_genes <- DE_OCvsRC
all_genes$SYMBOL <- all_genes$Geneid

# convert 
all_entrez <- bitr(all_genes$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
all_merged <- merge(all_entrez, all_genes, by.x = "SYMBOL", by.y = "Geneid")

# ranked vector by t stats
gene_list <- all_merged$t
names(gene_list) <- all_merged$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)  # rank from high to low

# remove duplicates 
sum(duplicated(names(gene_list))) 

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

# gseacurev1
gseaplot2(gsea_reactome, geneSetID = "R-HSA-9752946", title = "Expression and translocation of olfactory receptors")
# gseacurev 2
gseaplot2(gsea_reactome, geneSetID = "R-HSA-72695", title = "Formation of the ternary complex, and subsequently, the 43S complex")



################################

# OCCC vs GTEx normal ovary GSEA

all_ovary_genes <- DE_OCvsGTExOV
all_ovary_genes$SYMBOL <- all_ovary_genes$Geneid

# convert 
all_ovary_entrez <- bitr(all_ovary_genes$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
all_ovary_merged <- merge(all_ovary_entrez, all_ovary_genes, by.x = "SYMBOL", by.y = "Geneid")

# ranked vector by t
gene_list_ovary <- all_ovary_merged$t
names(gene_list_ovary) <- all_ovary_merged$ENTREZID
gene_list_ovary <- sort(gene_list_ovary, decreasing = TRUE)  # rank from high to low

# remove duplicates 
sum(duplicated(names(gene_list_ovary))) 
#gene_list_ovary <- gene_list_ovary[!duplicated(names(gene_list_ovary))]


# GO BP GSEA
gsea_go_ovary <- gseGO(geneList = gene_list_ovary,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

# KEGG GSEA
gsea_kegg_ovary <- gseKEGG(geneList = gene_list_ovary,
                     organism = "hsa",
                     minGSSize = 10,
                     pvalueCutoff = 0.05)

# Reactome GSEA
gsea_reactome_ovary <- gsePathway(geneList = gene_list_ovary,
                            organism = "human",
                            pvalueCutoff = 0.05,
                            verbose = FALSE)

dotplot(gsea_go_ovary, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_kegg_ovary, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_reactome_ovary, showCategory = 20, color = "NES") + theme_minimal()

# GSEA curves from reactome
#gseacurve3
gseaplot2(gsea_reactome_ovary, geneSetID = "R-HSA-373076", title = "Class A/1 (Rhodopsin-like receptors)")





###############################

###############################

# Kinase genes in OC vs RC

# kinase list ( threshold 0.05)
kinases <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/kinase_basic.xlsx", col_names = TRUE)
colnames(kinases)

kinases$Offical_gene_symbol <- toupper(kinases$Offical_gene_symbol)

# significant genes -- kinase list

up_kinases <- up_genes %>% filter(SYMBOL %in% kinases$Offical_gene_symbol)
down_kinases <- down_genes %>% filter(SYMBOL %in% kinases$Offical_gene_symbol)


up_kinases_entrez <- bitr(up_kinases$SYMBOL,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

down_kinases_entrez <- bitr(down_kinases$SYMBOL,
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = org.Hs.eg.db)

up_kinases_ids <- unique(up_kinases_entrez$ENTREZID)
down_kinases_ids <- unique(down_kinases_entrez$ENTREZID)

kegg_up_kinase <- enrichKEGG(gene = up_kinases_ids,
                             organism = "hsa",
                             pvalueCutoff = 0.05)

kegg_down_kinase <- enrichKEGG(gene = down_kinases_ids,
                               organism = "hsa",
                               pvalueCutoff = 0.05)
###
reactome_up_kinase <- enrichPathway(gene = up_kinases_ids,
                                    organism = "human",
                                    pvalueCutoff = 0.05)

reactome_down_kinase <- enrichPathway(gene = down_kinases_ids,
                                      organism = "human",
                                      pvalueCutoff = 0.05)
###

barplot(kegg_up_kinase, showCategory = 20) + theme_minimal()
barplot(reactome_up_kinase, showCategory = 20) + theme_minimal()

barplot(kegg_down_kinase, showCategory = 20) + theme_minimal()
barplot(reactome_down_kinase, showCategory = 20)+ theme_minimal()






