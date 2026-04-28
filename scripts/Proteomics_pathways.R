

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
#MSigDB
hallmark_gmt <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/h.all.v2026.1.Hs.symbols.gmt")
c7_msig <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/c7.all.v2026.1.Hs.symbols.gmt") # immunological
c6_msig <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/c6.all.v2026.1.Hs.symbols.gmt") # oncogenic
c3_msig <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/c3.tft.v2026.1.Hs.symbols.gmt") # TFT

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

#################
# msigdb
hallmark <- data.frame(gs_name = hallmark_gmt$term, gene_symbol = hallmark_gmt$gene)

all_merged_unique <- results_mapped %>%
  group_by(Geneid) %>%
  slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list_symbol <- all_merged_unique$logFC
names(gene_list_symbol) <- all_merged_unique$Geneid
gene_list_symbol <- sort(gene_list_symbol, decreasing = TRUE)
sum(duplicated(names(gene_list_symbol)))  


# hallmark
gsea_hallmark <- GSEA(geneList = gene_list_symbol,
                      TERM2GENE = hallmark,
                      minGSSize = 10,
                      pvalueCutoff = 0.05,
                      verbose = FALSE)

dotplot(gsea_hallmark, showCategory = 20, color = "NES") + theme_minimal()
cnetplot(gsea_hallmark, foldChange = gene_list_symbol, showCategory = 6)


# barplot
hallmark_results <- as.data.frame(gsea_hallmark) %>% mutate(
  direction = ifelse(NES > 0, "Up", "Down"),
  significance = case_when(p.adjust < 0.001 ~ "***",
                           p.adjust < 0.01 ~ "**",
                           p.adjust < 0.05 ~ "*",
                           TRUE ~ ""),
  Description = gsub("HALLMARK_", "", Description),
  Description = gsub("_", " ", Description)) %>% arrange(NES)  

ggplot(hallmark_results, aes(x = 1, y = reorder(Description, NES), fill = NES)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = significance), size = 5, vjust = 0.75) +
  scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "firebrick",
    midpoint = 0, name = "NES"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title   = element_blank(),
    panel.grid   = element_blank())




# save gsea
write_xlsx(list(GO_BP = as.data.frame(gsea_go),
                GO_MF = as.data.frame(gsea_mf),
                KEGG = as.data.frame(gsea_kegg),
                Reactome = as.data.frame(gsea_reactome),
                Hallmark_MSigDB = as.data.frame(gsea_hallmark)), "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/GSEA_Prot_OCvsRC_results.xlsx")




##########################################
# c6_msig - oncogenic pathway
##########################################

gsea_c6 <- GSEA(
  geneList  = gene_list_symbol,
  TERM2GENE = c6_msig,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff= 0.05,
  verbose = FALSE,
  seed= TRUE
)

head(as.data.frame(gsea_c6))

dotplot(gsea_c6, showCategory = 20, color = "NES") + theme_minimal() 
dotplot(gsea_c6, split = ".sign") + facet_grid(.~.sign)

cnetplot(gsea_c6, foldChange = gene_list_symbol, showCategory = 5)


############################################
# OCCC vs ccRCC c7_msig - immunological pathways/signatures
#############################################

gsea_c7 <- GSEA(
  geneList  = gene_list_symbol,
  TERM2GENE = c7_msig,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff= 0.05,
  verbose = FALSE,
  seed= TRUE
)

head(as.data.frame(gsea_c7))
cnetplot(gsea_c7, foldChange = gene_list_symbol, showCategory = 5)
dotplot(gsea_c7, split = ".sign") + facet_grid(.~.sign)

dotplot(gsea_c7, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_c7, showCategory = 20, color = "NES", split = ".sign") +
  facet_grid(. ~ .sign) +
  theme_minimal()

gseaplot2(gsea_c7, geneSetID = 1:5)


#################
# OCCC vs ccRCC - C3 TFT 
#################
gsea_c3 <- GSEA(
  geneList= gene_list_symbol,
  TERM2GENE = c3_msig,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff= 0.5,
  verbose = FALSE,
  seed= TRUE
)
# appears only at cut of value of 0.5 there is none significant 
head(as.data.frame(gsea_c3)) # none
dotplot(gsea_c3, split = ".sign") + facet_grid(.~.sign)
dotplot(gsea_c3, showCategory = 20, color = "NES") + theme_minimal()

# CEBP
######################

# KINASE

######################
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











