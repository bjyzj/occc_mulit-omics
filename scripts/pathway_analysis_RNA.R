
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(readxl)
library(writexl)
library(msigdbr)
suppressPackageStartupMessages(library(ExperimentHub))
suppressPackageStartupMessages(library(GSEABase))

set.seed(123)
# Load the data
DE_OCvsRC <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/transcriptomics_results/DE_results_OCCC_vs_ccRCC_ratio.csv")
DE_OCvsGTExOV <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/transcriptomics_results/DE_results_OCCC_vs_GTEx_Ovary_ratio.csv")
#MSigDB
hallmark <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/h.all.v2026.1.Hs.symbols.gmt")
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

#save
write_xlsx(list(KEGG_up = as.data.frame(kegg_up),
                KEGG_down = as.data.frame(kegg_down),
                Reactome_up = as.data.frame(reactome_up),
                Reactome_down = as.data.frame(reactome_down)), "ORA_results.xlsx")
                                
####################
# gene ontology ORA
ego_up <- enrichGO(gene = up_entrez_ids,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pvalueCutoff = 0.05)

ego_down <- enrichGO(gene = down_entrez_ids,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pvalueCutoff = 0.05)

# up
p_up <- barplot(ego_up,x = "Count",
                showCategory = 15,
                split = "ONTOLOGY"
) +
  aes(fill = ONTOLOGY) +
  scale_fill_brewer(palette = "Set2") +
  autofacet()

p_up + geom_text(aes(label = Count), hjust = -0.2, size = 3) +
  scale_y_discrete() +
  theme_classic() 

# down
p_down <- barplot(ego_down,x = "Count",
                  showCategory = 15,
                  split = "ONTOLOGY"
) +
  aes(fill = ONTOLOGY) +
  scale_fill_brewer(palette = "Set2") +
  autofacet()

p_down + geom_text(aes(label = Count), hjust = -0.2, size = 3) +
  scale_y_discrete() +
  theme_classic()

#####################

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
# GO MF
gsea_mf <- gseGO(geneList = gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont = "MF",
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
dotplot(gsea_mf, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_kegg, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_reactome, showCategory = 20, color = "NES") + theme_minimal()

# gseacurev1
gseaplot2(gsea_reactome, geneSetID = "R-HSA-9752946", title = "Expression and translocation of olfactory receptors")
# gseacurev 2
gseaplot2(gsea_reactome, geneSetID = "R-HSA-72695", title = "Formation of the ternary complex, and subsequently, the 43S complex")

gseaplot2(gsea_reactome, geneSetID = 1:5)


# msigdb gsea
# SYMBOL-based gene list (for Hallmark)
# symbol duplicated 2
all_merged_unique <- all_merged %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list_symbol <- all_merged_unique$t
names(gene_list_symbol) <- all_merged_unique$SYMBOL
gene_list_symbol <- sort(gene_list_symbol, decreasing = TRUE)
sum(duplicated(names(gene_list_symbol)))
# hallmark
gsea_hallmark <- GSEA(geneList = gene_list_symbol,
                      TERM2GENE = hallmark,
                      minGSSize = 10,
                      pvalueCutoff = 0.05,
                      verbose = FALSE)

dotplot(gsea_hallmark, showCategory = 20, color = "NES") + theme_minimal()
cnetplot(gsea_hallmark, foldChange = gene_list_symbol, showCategory = 3)
heatplot(gsea_hallmark, foldChange = gene_list_symbol, showCategory = 5)
heatplot(gsea_hallmark,
         foldChange = gene_list_symbol,
         showCategory = c(
           "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           "HALLMARK_TGF_BETA_SIGNALING",
           "HALLMARK_KRAS_SIGNALING_DN"
         ))


###############
# heatmap directionality on hallmark paths

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
  geom_text(aes(label = significance), size = 4, vjust = 0.75) +
  scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "firebrick",
    midpoint = 0, name = "NES"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title   = element_blank(),
    panel.grid   = element_blank()
  ) +
  labs(title = "Hallmark Pathways – OC vs RC")


# save gsea
write_xlsx(list(GO_BP = as.data.frame(gsea_go),
                GO_MF = as.data.frame(gsea_mf),
                KEGG = as.data.frame(gsea_kegg),
                Reactome = as.data.frame(gsea_reactome),
                Hallmark_MSigDB = as.data.frame(gsea_hallmark)), "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/GSEA_RNA_results_OCvsRC.xlsx")

###########
# gene + patwhays relationship heatmap
###########
gsea_go_readable <- setReadable(
  gsea_go,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

heatplot(gsea_go_readable, foldChange = gene_list, showCategory = 5)
cnetplot(gsea_go_readable, foldChange = gene_list, showCategory = 5)

gsea_kegg_readable <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
gsea_reactome_readable <- setReadable(gsea_reactome, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

heatplot(gsea_kegg_readable, foldChange = gene_list, showCategory = 5)
cnetplot(gsea_kegg_readable, foldChange = gene_list, showCategory = 5)

heatplot(gsea_reactome_readable, foldChange = gene_list, showCategory = 5)
cnetplot(gsea_reactome_readable, foldChange = gene_list, showCategory = 5)


# extracting common pathways of OCCC and ccRCC (NES=0)

go_res <- as.data.frame(gsea_go)
kegg_res <- as.data.frame(gsea_kegg)
reactome_res <- as.data.frame(gsea_reactome)

shared_go <- go_res %>%
  filter(p.adjust < 0.05 & abs(NES) <=0.5)

shared_kegg <- kegg_res %>%
  filter(p.adjust < 0.05 & abs(NES) < 1)

shared_reactome <- reactome_res %>%
  filter(p.adjust < 0.05 & abs(NES) < 1)


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

gsea_mf_ovary <- gseGO(geneList = gene_list_ovary,
                 OrgDb = org.Hs.eg.db,
                 ont = "MF",
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
dotplot(gsea_mf_ovary, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_kegg_ovary, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_reactome_ovary, showCategory = 20, color = "NES") + theme_minimal()

# GSEA curves from reactome
#gseacurve3
gseaplot2(gsea_reactome_ovary, geneSetID = "R-HSA-373076", title = "Class A/1 (Rhodopsin-like receptors)")


# msigdb gsea ovary

all_ovary_merged_unique <- all_ovary_merged %>%
  group_by(SYMBOL) %>%
  slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  ungroup()

gene_list_ovary_symbol <- all_ovary_merged_unique$t
names(gene_list_ovary_symbol) <- all_ovary_merged_unique$SYMBOL
gene_list_ovary_symbol <- sort(gene_list_ovary_symbol, decreasing = TRUE)
sum(duplicated(names(gene_list_ovary_symbol))) 

gsea_hallmark_ovary <- GSEA(geneList = gene_list_ovary_symbol,
                            TERM2GENE = hallmark,
                            minGSSize = 10,
                            pvalueCutoff = 0.05,
                            verbose = FALSE)

# heatmap
hallmark_ovary_results <- as.data.frame(gsea_hallmark_ovary) %>%
  mutate(direction = ifelse(NES > 0, "Up", "Down"),
         significance = case_when(p.adjust < 0.001 ~ "***",
                                  p.adjust < 0.01 ~ "**",
                                  p.adjust < 0.05 ~ "*",
                                  TRUE ~ ""),
         Description = gsub("HALLMARK_", "", Description),
         Description = gsub("_", " ", Description)) %>% arrange(NES)

ggplot(hallmark_ovary_results, aes(x = 1, y = reorder(Description, NES), fill = NES)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = significance), size = 4, vjust = 0.75) +
  scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "firebrick",
    midpoint = 0, name = "NES"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title   = element_blank(),
    panel.grid   = element_blank()
  ) +
  labs(title = "Hallmark Pathways – OCCC vs Normal Ovary")



write_xlsx(list(GO_BP = as.data.frame(gsea_go_ovary),
                GO_MF = as.data.frame(gsea_mf_ovary),
                KEGG = as.data.frame(gsea_kegg_ovary),
                Reactome = as.data.frame(gsea_reactome_ovary),
                Hallmark_MSigDB = as.data.frame(gsea_hallmark_ovary)), "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/GSEA_OCvsOVARY_RNA_results.xlsx")

###########################

# for shared pathways

go_ccrcc   <- as.data.frame(gsea_go)
go_ovary  <- as.data.frame(gsea_go_ovary)

kegg_ccrcc <- as.data.frame(gsea_kegg)
kegg_ovary <- as.data.frame(gsea_kegg_ovary)

reactome_ccrcc <- as.data.frame(gsea_reactome)
reactome_ovary  <- as.data.frame(gsea_reactome_ovary)

# shared significant pathways 
shared_go <- intersect(
  go_ccrcc$Description[go_ccrcc$p.adjust < 0.05],
  go_ovary$Description[go_ovary$p.adjust < 0.05]
)

shared_kegg <- intersect(
  kegg_ccrcc$Description[kegg_ccrcc$p.adjust < 0.05],
  kegg_ovary$Description[kegg_ovary$p.adjust < 0.05]
)

shared_reactome <- intersect(
  reactome_ccrcc$Description[reactome_ccrcc$p.adjust < 0.05],
  reactome_ovary$Description[reactome_ovary$p.adjust < 0.05]
)

# Print shared pathways
shared_go
shared_kegg
shared_reactome


# shared kegg
gsea_kegg_shared <- gsea_kegg[gsea_kegg$Description %in% shared_kegg, ]
# Keep only the pathways in shared_kegg
shared_idx <- which(gsea_kegg@result$Description %in% shared_kegg)
gsea_kegg_shared_obj <- gsea_kegg
gsea_kegg_shared_obj@result <- gsea_kegg@result[shared_idx, ]

# Make it readable for plotting
gsea_kegg_shared_readable <- setReadable(
  gsea_kegg_shared_obj,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

heatplot(gsea_kegg_shared_readable, foldChange = gene_list, showCategory = 10)
cnetplot(gsea_kegg_shared_readable, foldChange = gene_list, showCategory = 10)

# shared reactome
shared_idx_reactome <- which(gsea_reactome@result$Description %in% shared_reactome)
gsea_reactome_shared_obj <- gsea_reactome
gsea_reactome_shared_obj@result <- gsea_reactome@result[shared_idx_reactome, ]

gsea_reactome_shared_readable <- setReadable(
  gsea_reactome_shared_obj,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

heatplot(gsea_reactome_shared_readable, foldChange = gene_list, showCategory = 5)
cnetplot(gsea_reactome_shared_readable, foldChange = gene_list, showCategory = 10)



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


write_xlsx(list(Up_kinases = up_kinases,
                Down_kinases = down_kinases,
                KEGG_up_kinase = as.data.frame(kegg_up_kinase),
                KEGG_down_kinase = as.data.frame(kegg_down_kinase),
                Reactome_up_kinase = as.data.frame(reactome_up_kinase),
                Reactome_down_kinase = as.data.frame(reactome_down_kinase)), "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/Kinase_analysis.xlsx")



