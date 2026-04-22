
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)
library(readxl)
library(writexl)
library(pheatmap)
library(circlize)
#library(msigdbr)
#suppressPackageStartupMessages(library(ExperimentHub))
#suppressPackageStartupMessages(library(GSEABase))

set.seed(123)
# Load the data
DE_OCvsRC <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/transcriptomics_results/DE_results_OCCC_vs_ccRCC_ratio.csv")
DE_OCvsGTExOV <- read_csv("/Users/beyzaerkal/Desktop/occc_multi-omics/results/transcriptomics_results/DE_results_OCCC_vs_GTEx_Ovary_ratio.csv")
#MSigDB
hallmark <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/h.all.v2026.1.Hs.symbols.gmt")
c7_msig <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/c7.all.v2026.1.Hs.symbols.gmt") # immunological
c6_msig <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/c6.all.v2026.1.Hs.symbols.gmt") # oncogenic
c3_msig <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/c3.tft.v2026.1.Hs.symbols.gmt") # TFT

#####################
# OC vs RC ORA
ggplot(DE_OCvsGTExOV, aes(x = logFC)) +
  geom_histogram(bins = 50) +
  theme_minimal() +
  labs(x = "logFC", y = "Gene count", title = "Distribution of logFC (OCCC vs GTEx ovary)")

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
                Reactome_down = as.data.frame(reactome_down)), "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/ORA_RNA_OCvsRC_results.xlsx")
                                
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

dotplot(gsea_kegg, showCategory = 10, color = "NES") + theme_minimal()
dotplot(gsea_reactome, showCategory = 10, color = "NES") + theme_minimal()

# gseacurev1
gseaplot2(gsea_reactome, geneSetID = "R-HSA-9752946", title = "Expression and translocation of olfactory receptors")
# gseacurev 2
gseaplot2(gsea_reactome, geneSetID = "R-HSA-72695", title = "Formation of the ternary complex, and subsequently, the 43S complex")

gseaplot2(gsea_reactome, geneSetID = 1:5)

#############
# simplify
gsea_go_simplified <- simplify(gsea_go, cutoff = 0.7, by = "p.adjust", select_fun = min)
gsea_mf_simplified <- simplify(gsea_mf, cutoff = 0.7, by = "p.adjust", select_fun = min)

dotplot(gsea_go_simplified, showCategory = 20, color = "NES") + theme_minimal()
dotplot(gsea_mf_simplified, showCategory = 20, color = "NES") + theme_minimal()

gsea_go_simplified <- pairwise_termsim(gsea_go_simplified)
emapplot(gsea_go_simplified, showCategory = 30)

emapplot(pairwise_termsim(gsea_reactome), showCategory=30)

########
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
cnetplot(gsea_hallmark, foldChange = gene_list_symbol, showCategory = 6)
heatplot(gsea_hallmark, foldChange = gene_list_symbol, showCategory = 5)
heatplot(gsea_hallmark,
         foldChange = gene_list_symbol,
         showCategory = c(
           "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           "HALLMARK_TGF_BETA_SIGNALING",
           "HALLMARK_KRAS_SIGNALING_DN"
         ))
# shorter
top_genes <- sort(abs(gene_list_symbol), decreasing = TRUE)[1:100]
heatplot(gsea_hallmark,
         foldChange = gene_list_symbol[names(top_genes)],
         showCategory = c(
           "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
           "HALLMARK_TGF_BETA_SIGNALING",
           "HALLMARK_KRAS_SIGNALING_DN"
         ))

summary(gene_list_symbol[names(top_genes)])

# gene info 
gseaplot2(gsea_hallmark, geneSetID = 1:5)
gseaplot2(gsea_hallmark, geneSetID = 1:3)

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


#############
# leading edge genes - core enrichemnt
kegg_res <- as.data.frame(gsea_kegg)
kegg_res %>% dplyr::select(ID, Description, NES, p.adjust, core_enrichment)
# one path
ids <- strsplit(kegg_res$core_enrichment[1], "/")[[1]]

bitr(ids,
     fromType = "ENTREZID",
     toType = "SYMBOL",
     OrgDb = org.Hs.eg.db)

# get top negative 
neg_kegg <- kegg_res %>%
  filter(NES < 0, p.adjust < 0.05)

neg_kegg %>% dplyr::select(Description, NES, core_enrichment)

neg_genes <- neg_kegg %>%
  as_tibble() %>%
  dplyr::select(Description, core_enrichment) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(ENTREZID = core_enrichment)


symbols <- bitr(neg_genes$ENTREZID,
                fromType="ENTREZID",
                toType="SYMBOL",
                OrgDb=org.Hs.eg.db)

neg_genes <- left_join(neg_genes, symbols, by="ENTREZID")

# recurrign drivers
neg_genes %>%
  count(SYMBOL, sort = TRUE)

as.data.frame(gsea_kegg) %>%
  filter(NES < 0, p.adjust < 0.05) %>%
  select(Description, NES, core_enrichment)

neg_genes %>%
  count(SYMBOL, sort=TRUE) %>%
  slice_head(n=20) %>%
  ggplot(aes(n, reorder(SYMBOL, n))) +
  geom_col(fill="steelblue") +
  labs(x="Pathway recurrence", y="Gene")


# check the genes that appeared if tehya re diffenrtilly expressed
DE_OCvsRC %>%
  filter(Geneid %in% c("PRKACA","AKT3","ADCY5"))
# candidate driver table
tribble(
  ~Gene, ~logFC, ~Role,
  "PRKACA",-1.42,"PKA catalytic signaling",
  "AKT3",-1.53,"PI3K/AKT signaling",
  "ADCY5",-1.01,"cAMP production"
)

# hallmark leadign edge
hallmark_res <- as.data.frame(gsea_hallmark)

neg_hallmark <- hallmark_res %>%
  filter(NES < 0, p.adjust < 0.05)

neg_hallmark_genes <- neg_hallmark %>%
  dplyr::select(Description, core_enrichment) %>%
  tidyr::separate_rows(core_enrichment, sep="/") %>%
  dplyr::rename(SYMBOL = core_enrichment)

# recurrign hallmakrs
neg_hallmark_genes %>%
  count(SYMBOL, sort=TRUE)

neg_hallmark_genes %>%
  count(SYMBOL, sort=TRUE) %>%
  slice_head(n=20) %>%
  ggplot(aes(n, reorder(SYMBOL, n))) +
  geom_col(fill="steelblue") +
  labs(x="Hallmark recurrence", y="Gene")

# reactoem leaidng edge


reactome_res <- as.data.frame(gsea_reactome)

neg_reactome <- reactome_res %>%
  filter(NES < 0, p.adjust < 0.05)

neg_reactome_genes <- neg_reactome %>%
  dplyr::select(Description, core_enrichment) %>%
  tidyr::separate_rows(core_enrichment, sep="/") %>%
  dplyr::rename(ENTREZID = core_enrichment)


symbols <- bitr(
  neg_reactome_genes$ENTREZID,
  fromType="ENTREZID",
  toType="SYMBOL",
  OrgDb=org.Hs.eg.db
)

neg_reactome_genes <- left_join(
  neg_reactome_genes,
  symbols,
  by="ENTREZID"
)

neg_reactome_genes %>%
  count(SYMBOL, sort=TRUE)

neg_reactome_genes %>%
  count(SYMBOL, sort=TRUE) %>%
  slice_head(n=20) %>%
  ggplot(aes(n, reorder(SYMBOL, n))) +
  geom_col(fill="steelblue") +
  labs(x="Hallmark recurrence", y="Gene")



####################
# GAINS - leading edge
####################


# kegg
# POSITIVE KEGG leading-edge genes (higher in OCCC vs ccRCC)

kegg_res <- as.data.frame(gsea_kegg)

kegg_res %>%
  dplyr::select(ID, Description, NES, p.adjust, core_enrichment)

# one positive pathway example
ids <- strsplit(
  kegg_res$core_enrichment[which(kegg_res$NES > 0)[1]],
  "/"
)[[1]]

bitr(ids,
     fromType = "ENTREZID",
     toType = "SYMBOL",
     OrgDb = org.Hs.eg.db)

pos_kegg <- kegg_res %>%
  dplyr::filter(NES > 0, p.adjust < 0.05)

pos_kegg %>%
  dplyr::select(Description, NES, core_enrichment)

# extract leading-edge genes
pos_genes <- pos_kegg %>%
  tibble::as_tibble() %>%
  dplyr::select(Description, core_enrichment) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(ENTREZID = core_enrichment)

symbols <- bitr(pos_genes$ENTREZID,
                fromType = "ENTREZID",
                toType = "SYMBOL",
                OrgDb = org.Hs.eg.db)

pos_genes <- dplyr::left_join(pos_genes, symbols, by = "ENTREZID")

# recurring positive drivers
pos_genes %>%
  count(SYMBOL, sort = TRUE)

as.data.frame(gsea_kegg) %>%
  dplyr::filter(NES > 0, p.adjust < 0.05) %>%
  dplyr::select(Description, NES, core_enrichment)

pos_genes %>%
  count(SYMBOL, sort = TRUE) %>%
  slice_head(n = 20) %>%
  ggplot(aes(n, reorder(SYMBOL, n))) +
  geom_col(fill = "firebrick") +
  labs(x = "Pathway recurrence", y = "Gene")

##############
#hallmakr positive
hallmark_res <- as.data.frame(gsea_hallmark)

pos_hallmark <- hallmark_res %>%
  dplyr::filter(NES > 0, p.adjust < 0.05)

pos_hallmark_genes <- pos_hallmark %>%
  dplyr::select(Description, core_enrichment) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(SYMBOL = core_enrichment)

pos_hallmark_genes %>%
  count(SYMBOL, sort = TRUE)

pos_hallmark_genes %>%
  count(SYMBOL, sort = TRUE) %>%
  slice_head(n = 20) %>%
  ggplot(aes(n, reorder(SYMBOL, n))) +
  geom_col(fill = "firebrick") +
  labs(x = "Hallmark recurrence (positive NES)", y = "Gene")

################
# reactome positive

reactome_res <- as.data.frame(gsea_reactome)

pos_reactome <- reactome_res %>%
  dplyr::filter(NES > 0, p.adjust < 0.05)

pos_reactome_genes <- pos_reactome %>%
  dplyr::select(Description, core_enrichment) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::rename(ENTREZID = core_enrichment)

symbols <- bitr(
  pos_reactome_genes$ENTREZID,
  fromType = "ENTREZID",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)

pos_reactome_genes <- left_join(
  pos_reactome_genes,
  symbols,
  by = "ENTREZID"
)

# recurring positive drivers
pos_reactome_genes %>%
  count(SYMBOL, sort = TRUE)

# top recurring genes plot
pos_reactome_genes %>%
  count(SYMBOL, sort = TRUE) %>%
  slice_head(n = 20) %>%
  ggplot(aes(n, reorder(SYMBOL, n))) +
  geom_col(fill = "firebrick") +
  labs(x = "Reactome recurrence (positive NES)", y = "Gene")

###############


###########

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


###########################
# gsea vertical col bar plots
###########################
go_df <- as.data.frame(gsea_go)[, c("ID", "Description", "NES", "p.adjust")]
go_df$DB <- "GO_BP"

kegg_df <- as.data.frame(gsea_kegg)[, c("ID", "Description", "NES", "p.adjust")]
kegg_df$DB <- "KEGG"

react_df <- as.data.frame(gsea_reactome)[, c("ID", "Description", "NES", "p.adjust")]
react_df$DB <- "REACTOME"

get_top_paths <- function(df, n = 5) {
  
  up <- df %>%
    filter(NES > 0) %>%
    arrange(p.adjust) %>%
    head(n)
  
  down <- df %>%
    filter(NES < 0) %>%
    arrange(p.adjust) %>%
    head(n)
  
  bind_rows(up, down)
}

top_go <- get_top_paths(go_df, 10)
top_kegg <- get_top_paths(kegg_df, 10)
top_reactome <- get_top_paths(react_df, 10)

top_go$DB <- "GO_BP"
top_kegg$DB <- "KEGG"
top_reactome$DB <- "REACTOME"

gsea_top_all <- bind_rows(top_go, top_kegg, top_reactome)

top_kegg$Direction <- ifelse(top_kegg$NES > 0,
                             "Upregulated in OCCC",
                             "Downregulated in OCCC")

ggplot(top_kegg,
       aes(x = NES,
           y = reorder(Description, NES),
           fill = Direction)) +
  geom_col(width = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + scale_fill_manual(values = c(
    "Upregulated in OCCC" = "firebrick",
    "Downregulated in OCCC" = "steelblue"
  )) +
  labs(x = "NES", y = "Pathway", fill = "Enrichment direction")

#reactome
top_reactome$Direction <- ifelse(top_reactome$NES > 0,
                             "Upregulated in OCCC",
                             "Downregulated in OCCC")

ggplot(top_reactome,
       aes(x = NES,
           y = reorder(Description, NES),
           fill = Direction)) +
  geom_col(width = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + scale_fill_manual(values = c(
    "Upregulated in OCCC" = "firebrick",
    "Downregulated in OCCC" = "steelblue"
  )) +
  labs(x = "NES", y = "Pathway", fill = "Enrichment direction")


top_go$Direction <- ifelse(top_go$NES > 0,
                                 "Upregulated in OCCC",
                                 "Downregulated in OCCC")

ggplot(top_go,
       aes(x = NES,
           y = reorder(Description, NES),
           fill = Direction)) +
  geom_col(width = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() + scale_fill_manual(values = c(
    "Upregulated in OCCC" = "firebrick",
    "Downregulated in OCCC" = "steelblue"
  )) +
  labs(x = "NES", y = "Pathway", fill = "Enrichment direction")

##########################################
# OCCC vs ccRCC c6_msig - oncogenic pathway
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

dotplot(gsea_c6, showCategory = 20, color = "NES") + theme_minimal() +
gseaplot2(gsea_c6, geneSetID = 1:5)

# cleaner patwhasy names
c6_res <- as.data.frame(gsea_c6) %>%
  mutate(
    Direction = ifelse(NES > 0, "Up in OCCC", "Down in OCCC"),
    Description = gsub("_", " ", Description)
  )

# top enriched 6 
top_c6 <- c6_res %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 15)

ggplot(top_c6,
       aes(NES, reorder(Description, NES), fill = Direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(y = "", x = "NES")

# core enriched genes
c6_res %>% dplyr::select(Description, NES, p.adjust, core_enrichment) %>% head()

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

dotplot(gsea_c7, showCategory = 20, color = "NES") + theme_minimal()
gseaplot2(gsea_c7, geneSetID = 1:5)

# cleaner pathwas
c7_res <- as.data.frame(gsea_c7) %>%
  mutate(
    Direction = ifelse(NES > 0, "Up in OCCC", "Down in OCCC"),
    Description = gsub("_", " ", Description)
  )

top_c7 <- c7_res %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 15)

ggplot(top_c7,
       aes(NES, reorder(Description, NES), fill = Direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(title = "Top Immunologic Signatures (C7)",
       y = "", x = "NES")

# core enriched genes
c7_res %>% dplyr::select(Description, NES, p.adjust, core_enrichment) %>% head()


#################
# OCCC vs ccRCC - C3 TFT 
#################
gsea_c3 <- GSEA(
  geneList= gene_list_symbol,
  TERM2GENE = c3_msig,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff= 0.05,
  verbose = FALSE,
  seed= TRUE
)

head(as.data.frame(gsea_c3))

dotplot(gsea_c3, showCategory = 20, color = "NES") + theme_minimal()
gseaplot2(gsea_c3, geneSetID = 1:5)

# cleaner pathwas
c3_res <- as.data.frame(gsea_c3) %>%
  mutate(
    Direction = ifelse(NES > 0, "Up in OCCC", "Down in OCCC"),
    Description = gsub("_", " ", Description)
  )

top_c3 <- c3_res %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 15)

ggplot(top_c3,
       aes(NES, reorder(Description, NES), fill = Direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(title = "Top Immunologic Signatures (C3)",
       y = "", x = "NES")

# core enriched genes
c3_res %>% dplyr::select(Description, NES, p.adjust, core_enrichment) %>% head()



write_xlsx(list(C6_Oncogenic_OCvsRC= as.data.frame(gsea_c6), C7_Immunologic_OCvsRC = as.data.frame(gsea_c7), C3_TFTs_OCvsRC = as.data.frame(gsea_c3)),
  "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/GSEA_C3_C6_C7_OCvsRC.xlsx")

################################

# OCCC vs GTEx normal ovary GSEA

################################
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

#########################
# msigdb gsea ovary
#########################

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

hallmark_ovary_results_plot <- hallmark_ovary_results %>%
  group_by(direction) %>%
  slice_max(order_by = abs(NES), n = 8) %>%
  ungroup()

ggplot(hallmark_ovary_results_plot, aes(x = 1, y = reorder(Description, NES), fill = NES)) +
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
    panel.grid   = element_blank())



write_xlsx(list(GO_BP = as.data.frame(gsea_go_ovary),
                GO_MF = as.data.frame(gsea_mf_ovary),
                KEGG = as.data.frame(gsea_kegg_ovary),
                Reactome = as.data.frame(gsea_reactome_ovary),
                Hallmark_MSigDB = as.data.frame(gsea_hallmark_ovary)), "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/GSEA_OCvsOVARY_RNA_results.xlsx")

########################

##########################################
# OCCC vs Ovary c6_msig - oncogenic pathway
##########################################

gsea_c6_ovary <- GSEA(
  geneList= gene_list_ovary_symbol,
  TERM2GENE = c6_msig,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff= 0.05,
  verbose = FALSE,
  seed= TRUE)

head(as.data.frame(gsea_c6_ovary))

dotplot(gsea_c6_ovary, showCategory = 20, color = "NES") + theme_minimal() 
gseaplot2(gsea_c6_ovary, geneSetID = 1:5)

# cleaner patwhasy names
c6_res_ovary <- as.data.frame(gsea_c6_ovary) %>%
  mutate(
    Direction = ifelse(NES > 0, "Up in OCCC", "Down in OCCC"),
    Description = gsub("_", " ", Description))

# top enriched 6 
top_c6_ovary <- c6_res_ovary %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 15)

ggplot(top_c6_ovary,
       aes(NES, reorder(Description, NES), fill = Direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(y = "", x = "NES")

# core enriched genes
c6_res_ovary %>% dplyr::select(Description, NES, p.adjust, core_enrichment) %>% head()

############################################
# OCCC vs Ovary c7_msig - immunological pathways/signatures
#############################################

gsea_c7_ovary <- GSEA(geneList = gene_list_ovary_symbol,
                TERM2GENE = c7_msig,
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff= 0.05,
                verbose = FALSE,
                seed= TRUE)

head(as.data.frame(gsea_c7_ovary))

dotplot(gsea_c7_ovary, showCategory = 20, color = "NES") + theme_minimal()
gseaplot2(gsea_c7_ovary, geneSetID = 1:5)

# cleaner pathwas
c7_res_ovary <- as.data.frame(gsea_c7_ovary) %>%
  mutate(Direction = ifelse(NES > 0, "Up in OCCC", "Down in OCCC"),
         Description = gsub("_", " ", Description))

top_c7_ovary <- c7_res_ovary %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 15)

ggplot(top_c7_ovary,
       aes(NES, reorder(Description, NES), fill = Direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(title = "Top Immunologic Signatures (C7)",
       y = "", x = "NES")

# core enriched genes
c7_res_ovary %>% dplyr::select(Description, NES, p.adjust, core_enrichment) %>% head()


#################
# OCCC vs Ovary - C3 TFT 
#################
gsea_c3_ovary <- GSEA(
  geneList= gene_list_ovary_symbol,
  TERM2GENE = c3_msig,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff= 0.05,
  verbose = FALSE,
  seed  = TRUE
)

head(as.data.frame(gsea_c3_ovary))

dotplot(gsea_c3_ovary, showCategory = 20, color = "NES") + theme_minimal()
gseaplot2(gsea_c3_ovary, geneSetID = 1:5)

# cleaner pathwas
c3_res_ovary <- as.data.frame(gsea_c3_ovary) %>%
  mutate(
    Direction = ifelse(NES > 0, "Up in OCCC", "Down in OCCC"),
    Description = gsub("_", " ", Description)
  )

top_c3_ovary <- c3_res_ovary %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 15)

ggplot(top_c3_ovary,
       aes(NES, reorder(Description, NES), fill = Direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = 2) +
  theme_bw() +
  labs(title = "Top Immunologic Signatures (C3)",
       y = "", x = "NES")

# core enriched genes
c3_res_ovary %>% dplyr::select(Description, NES, p.adjust, core_enrichment) %>% head()



write_xlsx(list(C6_Oncogenic_OCvsOvary= as.data.frame(gsea_c6_ovary), C7_Immunologic_OCvsOvary = as.data.frame(gsea_c7_ovary), 
                C3_TFTs_OCvsOvary = as.data.frame(gsea_c3_ovary)),
           "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/GSEA_C3_C6_C7_OCvsOvary.xlsx")



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

##############################

###############################

###############################

# Kinase genes in OC vs RC

###############################
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
                Reactome_down_kinase = as.data.frame(reactome_down_kinase)), "/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/Kinase_ORA_transcriptomics_analysis.xlsx")



#################
# circo plot with directionality
#################
up_df   <- as.data.frame(reactome_up_kinase)
down_df <- as.data.frame(reactome_down_kinase)

# top pathways 
up_df   <- up_df   %>% arrange(p.adjust) %>% dplyr::slice(1:10)
down_df <- down_df %>% arrange(p.adjust) %>% dplyr::slice(1:10)

# extract pathway- gene relationships 
extract_edges <- function(enrich_df, direction){
  
  all_edges <- data.frame()
  
  for(i in 1:nrow(enrich_df)){
    
    pathway_name <- enrich_df$Description[i]
    
    ids <- unlist(strsplit(enrich_df$geneID[i], "/"))
    
    conv <- bitr(ids,
                 fromType = "ENTREZID",
                 toType = "SYMBOL",
                 OrgDb = org.Hs.eg.db)
    
    temp <- data.frame(
      pathway   = pathway_name,
      kinase    = unique(conv$SYMBOL),
      direction = direction
    )
    
    all_edges <- rbind(all_edges, temp)
  }
  
  return(all_edges)
}

edges_up   <- extract_edges(up_df, "UP")
edges_down <- extract_edges(down_df, "DOWN")

edges <- bind_rows(edges_up, edges_down)

# only kinase genes
kinase_symbols <- unique(kinases$Offical_gene_symbol)

edges <- edges %>%
  filter(kinase %in% kinase_symbols)
# duplicates managed
edges <- unique(edges)

# most connected kinases (top 20)
top_kinases <- edges %>%
  count(kinase, sort = TRUE) %>%
  dplyr::slice(1:20) %>%
  pull(kinase)

edges <- edges %>%
  filter(kinase %in% top_kinases)
# colour
pathways <- unique(edges$pathway)
genes    <- unique(edges$kinase)

grid.col <- c(
  setNames(rep("grey40", length(pathways)), pathways),
  setNames(rep("firebrick", length(genes)), genes)
)


link_cols <- ifelse(edges$direction == "UP",
                    rgb(1,0,0,0.45),   # red
                    rgb(0,0,1,0.45))   # blue

circos.clear()

chordDiagram(
  x = edges[,c("pathway","kinase")],
  grid.col = grid.col,
  col = link_cols,
  annotationTrack = "grid",
  transparency = 0.25,
  preAllocateTracks = 1
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y){
    
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    circos.text(
      mean(xlim),
      ylim[1] + 0.1,
      sector.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0,0.5),
      cex = 0.55
    )
  },
  bg.border = NA
)



















