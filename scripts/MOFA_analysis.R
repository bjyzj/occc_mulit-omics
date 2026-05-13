
library(tidyverse)
library(ggplot2)
#library(ggrepel)
library(ComplexHeatmap)
#library(ggpubr) 
#library(matrixStats)
library(MOFA2)
#library(GGally)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(openxlsx)
library(cowplot)

# use docker mofa2/mofa2:latest
set.seed(23)

# MOFA only OCCC
RNA_seq_log2 <- readRDS("scripts/expression_matrix_occc_all_log2.rds")
Prot_log2 <- read_csv("scripts/proteomics_JJ2022_log2_centered.csv")

# curated pathway database
hallmark <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/h.all.v2026.1.Hs.symbols.gmt")
go_c5 <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/c5.all.v2026.1.Hs.symbols.gmt")

# metadata for RNA-seq
samples_RNA <- colnames(RNA_seq_log2)
ea_samples <- c("ES-2", "JHOC-5", "OVISE", "OVMANA", "OVTOKO", "RMG-I", "TOV-21G")

study_RNA <- ifelse(samples_RNA %in% ea_samples, "Expression Atlas (2026)",
                ifelse(samples_RNA %in% c("C1", "C2", "C3", "C4", "C5", "C6"), "Nagasawa et al. (2019)",
                       ifelse(grepl("^CCC_", samples_RNA), "GSE189553",
                              ifelse(grepl("^OVA[0-9]+", samples_RNA), "GSE160692",
                                     ifelse(grepl("IGO_07456", samples_RNA), "Bolton et al. (2022)", NA)))))


# metadata data.frame
sample_info_RNA <- data.frame(
  Sample = samples_RNA,
  Study = study_RNA,
  SourceType = ifelse(study_RNA == "Expression Atlas (2026)", "cell_line", "primary_tumour"),
  row.names = samples_RNA)



# metadata for Proteomics

samples_ProtOC <- colnames(Prot_log2)

sample_info_ProtOC <- data.frame(
  Sample = samples_ProtOC,
  Study = "Ji et al. (2022)",
  SourceType = "primary_tumour",
  row.names = samples_ProtOC
)

##########################
# combine all metadata
sample_info_all <- rbind(sample_info_RNA, sample_info_ProtOC)
# reorder to match combined sample order in views
all_samples <- c(colnames(RNA_seq_log2),
                 colnames(Prot_log2))
sample_info_all <- sample_info_all[all_samples, , drop = FALSE]

##############
# MOFA 
##############

# avoid full_join to introduce sourcetype like grouping but just geneid
# combine Geneid

Prot_log2 <- as.data.frame(Prot_log2)

colnames(Prot_log2)[1]  <- "Geneid"

sum(duplicated(Prot_log2$Geneid))
# check NA
sum(is.na(Prot_log2))

# set rownames and remove Geneid col
rownames(Prot_log2) <- Prot_log2$Geneid
Prot_log2 <- Prot_log2 %>% dplyr::select(-Geneid)
# matrix
Prot_log2 <- as.matrix(Prot_log2)

############
# create view list and align samples/paddings
###########

# make sure views match the metadata
# remove Geneid column and create list of views
views <- list(RNA = as.matrix(RNA_seq_log2),
              Proteomics = Prot_log2)
#rownames(views$RNA) <- RNA_seq_log2$Geneid

# checks alignment with metadata

head(RNA_seq_log2$Geneid)
dim(RNA_seq_log2)
rownames(RNA_seq_log2)[1:20]

# feature selection before MOFA object (recommended by MOFAobject message)
# align samples across views by padding with NAs - coz of differnet sample sets
all_samples <- unique(c(colnames(views$RNA), colnames(views$Proteomics)))

# Pad RNA
missing_RNA <- setdiff(all_samples, colnames(views$RNA))
if(length(missing_RNA) > 0){
  views$RNA <- cbind(
    views$RNA,
    matrix(NA, nrow = nrow(views$RNA), ncol = length(missing_RNA),
           dimnames = list(rownames(views$RNA), missing_RNA))
  )
}

# Pad Proteomics
missing_Prot <- setdiff(all_samples, colnames(views$Proteomics))
if(length(missing_Prot) > 0){
  views$Proteomics <- cbind(
    views$Proteomics,
    matrix(NA, nrow = nrow(views$Proteomics), ncol = length(missing_Prot),
           dimnames = list(rownames(views$Proteomics), missing_Prot))
  )
}

# Reorder columns to match
views$RNA <- views$RNA[, all_samples]
views$Proteomics <- views$Proteomics[, all_samples]

# Align metadata
sample_info_all <- sample_info_all[all_samples, , drop = FALSE]

all(colnames(views$RNA) %in% rownames(sample_info_all))
all(colnames(views$Proteomics) %in% rownames(sample_info_all))

setdiff(colnames(views$RNA), rownames(sample_info_all))
setdiff(colnames(views$Proteomics), rownames(sample_info_all))

# ------------- MOFA ---------------


num_fac<-c(10, 15, 20, 30) # number of factors to try - adjust based on data for variability

for(fac_idx in seq_along(num_fac)){
  
  # skipping ones that run before
  output_file <- file.path(getwd(),
                           paste("MOFA_output", as.character(num_fac[fac_idx]),"F.hdf5", sep=""))
  
  if (file.exists(output_file)) {
    message("Skipping training for ", num_fac[fac_idx], " factors as output file already exists.")
    next 
  }
  
  MOFAobject <- create_mofa(views, sample_metadata = sample_info_all) # data
  
  
  data_opts <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors<-num_fac[fac_idx]
  train_opts <- get_default_training_options(MOFAobject)
  
  MOFAobject <- prepare_mofa(object = MOFAobject,
                             data_options = data_opts,
                             model_options = model_opts,
                             training_options = train_opts
  )
  
  output_file = file.path(getwd(),
                          paste("MOFA_output", as.character(num_fac[fac_idx]),"F.hdf5", sep=""))
  # train
  MOFAobject_results <- run_mofa(MOFAobject, output_file, use_basilisk = TRUE)
  
  print(paste("Training complete for", num_fac, "factors."))
}


##########################
# MOFA results to visualise
m10 <- load_model("MOFA_results/MOFA_output10F.hdf5")
m15 <- load_model("MOFA_results/MOFA_output15F.hdf5")
m20 <- load_model("MOFA_results/MOFA_output20F.hdf5")
m30 <- load_model("MOFA_results/MOFA_output30F.hdf5")

slotNames(m20)
names(m20@data)
m20@covariates

get_elbo(m10)  # -1,430,339
get_elbo(m15)  # -575,056
get_elbo(m20)  # 51,892.38
get_elbo(m30)   # 866,758.1

# Total variance explained per view
head(get_variance_explained(m20)$r2_total[[1]])
# Variance explained for every factor in per view
head(get_variance_explained(m20)$r2_per_factor[[1]])

plot_variance_explained(m10)
plot_variance_explained(m15)
plot_variance_explained(m20) 
plot_variance_explained(m30) 


plot_data_overview(m20)
plot_data_overview(m30)

#samples_metadata(m20) |> colnames()
# adding more metadata labels
md <- samples_metadata(m20)

md_merged <- merge(md, 
                   sample_info_all,
                   by.x = "sample",
                   by.y = "Sample",
                   all.x = TRUE,
                   sort = FALSE)
samples_metadata(m20) <- md_merged
colnames(samples_metadata(m20))

plot_factor(m20, factors = 1:5, color_by = "Study", shape_by = "SourceType")

# Plot top weights
p1 <- plot_top_weights(m20, view = "RNA", factor = 1, nfeatures = 20, scale = TRUE, abs = FALSE) # not cancer specific
p2 <- plot_top_weights(m20, view = "RNA", factor = 2, nfeatures = 20, scale = TRUE, abs = FALSE) # signalling + immune rich epithelial
p3 <- plot_top_weights(m20, view = "RNA", factor = 3, nfeatures = 20, scale = TRUE, abs = FALSE) # mucin rich epithelial + secretory tumour phenotype
p4 <- plot_top_weights(m20, view = "RNA", factor = 4, nfeatures = 20, scale = TRUE, abs = FALSE) # tumour/stroma/immune interaction + fibrobalst + inflammation (microenv intergated interactions)

plot_grid(p1, p2, p3, p4, ncol = 2)

p5 <- plot_top_weights(m20, view = "Proteomics", factor = 1, nfeatures = 20, scale = TRUE, abs = FALSE) # cancer signatures
p6 <- plot_top_weights(m20, view = "Proteomics", factor = 2, nfeatures = 20, scale = TRUE, abs = FALSE) # non cancer specific: structural, rna processing
p7 <- plot_top_weights(m20, view = "Proteomics", factor = 3, nfeatures = 20, scale = TRUE, abs = FALSE) # hypoxia, metabolism, invasion related
p8 <- plot_top_weights(m20, view = "Proteomics", factor = 4, nfeatures = 20, scale = TRUE, abs = FALSE) # immune hot/poor, tumour stroma, ealry invasion (wnt), tumour microenv

plot_grid(p5, p6, p7, p8, ncol = 2)


# Violin plots for factors 1, 2, 3 coloured by Study
p_violin <- plot_factor(
  m20,
  factors = c(1, 2, 3, 4),
  color_by = "Study", 
  dot_size = 3,
  dodge = TRUE,
  legend = TRUE,
  add_violin = TRUE,
  violin_alpha = 0.25,
)
print(p_violin)

colnames(samples_metadata(m20))
plot_data_heatmap(m20, 
                  factor = 1, 
                  view = "RNA", 
                  features = 20,
                  denoise = TRUE,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = F, show_rownames = T,
                  annotation_samples = "Study",  
                  annotation_colors = list("Study"), 
                  annotation_legend = F,
                  scale = "row"
)

# checks: not a lot of correlation between factors meaning good model fit
plot_factor_cor(m20)


##############
# GSEA 
##############

# format hallmark- binary matrix version
genes <- unique(as.character(hallmark$gene))
terms <- unique(as.character(hallmark$term))

hallmark_mat <- matrix(0, nrow = length(terms), ncol = length(genes), dimnames = list(terms, genes))

for(i in seq_len(nrow(hallmark))) {
  hallmark_mat[
    as.character(hallmark$term[i]),
    as.character(hallmark$gene[i])
  ] <- 1
}

dim(hallmark_mat)

# c5
genes <- unique(as.character(go_c5$gene))
terms <- unique(as.character(go_c5$term))

go_c5_mat <- matrix(0, nrow = length(terms), ncol = length(genes), dimnames = list(terms, genes))

for(i in seq_len(nrow(go_c5))) {
  go_c5_mat[
    as.character(go_c5$term[i]),
    as.character(go_c5$gene[i])
  ] <- 1
}

dim(go_c5_mat)

##############
# GSEA RNA
##############

hallmark_pos <- run_enrichment(
  object = m20,
  view = "RNA",
  feature.sets = hallmark_mat,
  sign = "positive",
  statistical.test = "parametric"
)

hallmark_neg <- run_enrichment(
  object = m20,
  view = "RNA",
  feature.sets = hallmark_mat,
  sign = "negative",
  statistical.test = "parametric"
)

hallmark_pos$sigPathways[[1]]
hallmark_neg$sigPathways[[1]]
names(hallmark_pos)
hallmark_pos$set.statistics[1:5,1]
hallmark_neg$set.statistics[1:5,1]


# Factor 2
plot_enrichment(hallmark_pos, factor = 2,max.pathways = 15)
plot_enrichment(hallmark_neg, factor = 2,max.pathways = 15)

plot_enrichment_detailed(hallmark_pos, factor = 2,  max.genes = 8, max.pathways = 5)
plot_enrichment_detailed(hallmark_neg, factor = 2,  max.genes = 8, max.pathways = 5)
# Factor 3
plot_enrichment(hallmark_pos, factor = 3,max.pathways = 15)
plot_enrichment(hallmark_neg, factor = 3,max.pathways = 15)

plot_enrichment_detailed(hallmark_pos, factor = 3,  max.genes = 8, max.pathways = 5)
plot_enrichment_detailed(hallmark_neg, factor = 3,  max.genes = 8, max.pathways = 5)
# Factor 4
plot_enrichment(hallmark_pos, factor = 4,max.pathways = 15)
plot_enrichment(hallmark_neg, factor = 4,max.pathways = 15)

p_rna_4_pos <- plot_enrichment_detailed(hallmark_pos, factor = 4,  max.genes = 8, max.pathways = 5)
p_rna_4_neg <- plot_enrichment_detailed(hallmark_neg, factor = 4,  max.genes = 8, max.pathways = 5)



p_rna_4_pos +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )


p_rna_4_neg +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )


##################
# C5
c5_pos <- run_enrichment(
  object = m20,
  view = "RNA",
  feature.sets = go_c5_mat,
  sign = "positive",
  statistical.test = "parametric"
)

c5_neg <- run_enrichment(
  object = m20,
  view = "RNA",
  feature.sets = go_c5_mat,
  sign = "negative",
  statistical.test = "parametric"
)

c5_pos$sigPathways[[1]]
c5_neg$sigPathways[[1]]


# Factor 2
plot_enrichment(c5_pos, factor = 2,max.pathways = 15)
plot_enrichment(c5_neg, factor = 2,max.pathways = 15)

plot_enrichment(c5_pos, factor = 3,max.pathways = 15)
plot_enrichment(c5_neg, factor = 3,max.pathways = 15)

plot_enrichment(c5_pos, factor = 4,max.pathways = 15)
plot_enrichment(c5_neg, factor = 4,max.pathways = 15)




#############
# GSEA Prot
#############

hallmark_pos_prot <- run_enrichment(
  object = m20,
  view = "Proteomics",
  feature.sets = hallmark_mat,
  sign = "positive",
  statistical.test = "parametric"
)

hallmark_neg_prot <- run_enrichment(
  object = m20,
  view = "Proteomics",
  feature.sets = hallmark_mat,
  sign = "negative",
  statistical.test = "parametric"
)

hallmark_pos_prot$sigPathways[[1]]
hallmark_neg_prot$sigPathways[[1]]

hallmark_pos_prot$set.statistics[1:5,1]
hallmark_neg_prot$set.statistics[1:5,1]

# none for factor 1,2, 5,6

plot_enrichment(hallmark_pos_prot, factor = 3,max.pathways = 15)
plot_enrichment_detailed(hallmark_pos_prot, factor = 3,  max.genes = 8, max.pathways = 5)
plot_enrichment(hallmark_neg_prot, factor = 3,max.pathways = 15) # none

plot_enrichment(hallmark_pos_prot, factor = 4,max.pathways = 15)
p_prot_f4 <- plot_enrichment_detailed(hallmark_pos_prot, factor = 4,  max.genes = 8, max.pathways = 5)
plot_enrichment(hallmark_neg_prot, factor = 4,max.pathways = 15) # none
 
p_prot_f4 + theme_bw(base_size = 13) +
  theme(
    axis.text.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )


# C5
c5_pos_prot <- run_enrichment(
  object = m20,
  view = "Proteomics",
  feature.sets = go_c5_mat,
  sign = "positive",
  statistical.test = "parametric"
)

c5_neg_prot <- run_enrichment(
  object = m20,
  view = "Proteomics",
  feature.sets = go_c5_mat,
  sign = "negative",
  statistical.test = "parametric"
)

c5_pos_prot$sigPathways[[1]]
c5_neg_prot$sigPathways[[1]]

c5_pos_prot$set.statistics[1:5,1]
c5_neg_prot$set.statistics[1:5,1]

plot_enrichment_detailed(c5_pos_prot, factor = 3)
plot_enrichment_detailed(c5_neg_prot, factor = 3,  max.genes = 8, max.pathways = 5)


## save results
########################

rna_factor1_df <- get_weights(m20, view = "RNA", factor = 1, as.data.frame = TRUE)
rna_factor2_df <- get_weights(m20, view = "RNA", factor = 2, as.data.frame = TRUE)
rna_factor3_df <- get_weights(m20, view = "RNA", factor = 3, as.data.frame = TRUE)
rna_factor4_df <- get_weights(m20, view = "RNA", factor = 4, as.data.frame = TRUE)

prot_factor1_df <- get_weights(m20, view = "Proteomics", factor = 1, as.data.frame = TRUE)
prot_factor2_df <- get_weights(m20, view = "Proteomics", factor = 2, as.data.frame = TRUE)
prot_factor3_df <- get_weights(m20, view = "Proteomics", factor = 3, as.data.frame = TRUE)
prot_factor4_df <- get_weights(m20, view = "Proteomics", factor = 4, as.data.frame = TRUE)


wb <- createWorkbook()

# RNA sheets
addWorksheet(wb, "RNA_F1")
writeData(wb, "RNA_F1", rna_factor1_df)

addWorksheet(wb, "RNA_F2")
writeData(wb, "RNA_F2", rna_factor2_df)

addWorksheet(wb, "RNA_F3")
writeData(wb, "RNA_F3", rna_factor3_df)

addWorksheet(wb, "RNA_F4")
writeData(wb, "RNA_F4", rna_factor4_df)


# Proteomics sheets
addWorksheet(wb, "Prot_F1")
writeData(wb, "Prot_F1", prot_factor1_df)

addWorksheet(wb, "Prot_F2")
writeData(wb, "Prot_F2", prot_factor2_df)

addWorksheet(wb, "Prot_F3")
writeData(wb, "Prot_F3", prot_factor3_df)

addWorksheet(wb, "Prot_F4")
writeData(wb, "Prot_F4", prot_factor4_df)

saveWorkbook(wb,"Supplementary_MOFA_weights.xlsx",overwrite = TRUE)








