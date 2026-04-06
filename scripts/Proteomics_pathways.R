

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
library(writexl)
set.seed(123)


## Load data
DEProt <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/proteomics_results/DE_OCCC_vs_ccRCC_proteomics_qn.csv")

# ORA
results <- DEProt 

gene_map <- bitr(results$Geneid,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

results_mapped <- results %>% inner_join(gene_map, by = c("Geneid" = "SYMBOL"))

up_genes <- results_mapped %>% filter(adj.P.Val < 0.05 & logFC > 0.25)
down_genes <- results_mapped %>% filter(adj.P.Val < 0.05 & logFC < -0.25)

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


##########


sig_genes <- results_mapped %>%
  filter(adj.P.Val < 0.05) %>%
  pull(ENTREZID) %>%
  unique()

ego_all <- enrichGO(
  gene = sig_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",             
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

p_nofacet <- barplot(
  ego_all,
  x = "Count",
  showCategory = 15,
  split = "ONTOLOGY"
) +
  aes(fill = ONTOLOGY) +
  scale_fill_brewer(palette = "Set2")

p_nofacet +
  geom_text(aes(label = Count), hjust = -0.2, size = 3) +
  theme_classic() +
  labs(x = "Gene Count", y = NULL)






########


# GSEA

gene_list <- results_mapped$t   
names(gene_list) <- results_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

# GO BP
gsea_go <- gseGO(geneList = gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont = "BP",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

gsea_mf <- gseGO(geneList = gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont = "MF",
                 minGSSize = 10,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05,
                 verbose = FALSE)

gsea_kegg <- gseKEGG(geneList = gene_list, 
                     organism = "hsa", 
                     minGSSize = 10, 
                     pvalueCutoff = 0.05)

gsea_reactome <- gsePathway(geneList = gene_list,
                            organism = "human", 
                            pvalueCutoff = 0.05, verbose = FALSE)


dotplot(gsea_go, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_mf, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_kegg, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_reactome, showCategory = 20, color = "NES") + theme_minimal()

#require(DOSE)
dotplot(gsea_reactome, showCategory=10, split=".sign") + facet_grid(.~.sign)

# save gsea
write_xlsx(list(GO_BP = as.data.frame(gsea_go),
                GO_MF = as.data.frame(gsea_mf),
                KEGG = as.data.frame(gsea_kegg),
                Reactome = as.data.frame(gsea_reactome)), "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/GSEA_Prot_OCvsRC_results.xlsx")



# KINASE


# kinase list ( threshold 0.05)
kinases <- read_excel("/Users/beyzaerkal/Desktop/internship/internship_env/kinase_basic.xlsx", col_names = TRUE)
colnames(kinases)

# SYMBOL col match kinase list
up_genes$SYMBOL <- toupper(up_genes$Geneid)
down_genes$SYMBOL <- toupper(down_genes$Geneid)

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











