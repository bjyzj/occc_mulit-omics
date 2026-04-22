
library(tidyverse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(GSEABase)
library(writexl)
library(readxl)
library(circlize)

# hallmark 
hallmark <- read.gmt("/Users/beyzaerkal/Desktop/internship/internship_env/h.all.v2026.1.Hs.symbols.gmt")

results_ssgsea <- read.csv("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/ssgsea_scores.csv", row.names = 1)
# output from ecotyper

process_ecotyper_nested <- function(root_path) {
  
  cell_types <- c(
    "B.cells",
    "CD4.T.cells",
    "CD8.T.cells",
    "Dendritic.cells",
    "Endothelial.cells",
    "Epithelial.cells",
    "Fibroblasts",
    "Mast.cells",
    "Monocytes.and.Macrophages",
    "NK.cells",
    "PCs",
    "PMNs"
  )
  
  samples <- list.files(root_path, full.names = TRUE)
  all_results <- list()
  
  for (sample_path in samples) {
    
    cell_state_path <- file.path(sample_path, "Carcinoma_Cell_States")
    if (!dir.exists(cell_state_path)) next
    
    for (ct in cell_types) {
      
      ct_path <- file.path(cell_state_path, ct)
      if (!dir.exists(ct_path)) next
      
      abd_file<- list.files(ct_path, pattern = "Abundance",  full.names = TRUE)
      assn_file <- list.files(ct_path, pattern = "Assignment", full.names = TRUE)
      
      if (length(abd_file) == 0 || length(assn_file) == 0) next
      
      abd <- tryCatch(
        read.delim(abd_file[1],  check.names = FALSE),
        error = function(e) NULL
      )
      
      assn <- tryCatch(
        read.delim(assn_file[1], check.names = FALSE),
        error = function(e) NULL
      )
      
      if (is.null(abd) || is.null(assn)) next
      if (!is.data.frame(abd) || !is.data.frame(assn)) next
      if (ncol(abd) < 2 || ncol(assn) < 2) next
      
      colnames(assn) <- trimws(gsub("\\s+", " ", colnames(assn)))
      colnames(abd) <- trimws(gsub("\\s+", " ", colnames(abd)))
      
      colnames(abd)[1]<- "ID"
      colnames(assn)[1] <- "ID"
      
      if (!any(grepl("Cell", colnames(assn)))) next
      cell_col <- grep("Cell", colnames(assn), value = TRUE)[1]
      
      assn_small <- assn[, c("ID", cell_col)]
      colnames(assn_small)[2] <- "Cell.State"
      
      abd_long <- melt(
        abd,
        id.vars = "ID",
        variable.name = "State",
        value.name = "Abundance"
      )
      
      df_long <- merge(abd_long, assn_small, by = "ID", all.x = TRUE)
      
      df_long$CellType <- ct
      df_long$Sample <- basename(sample_path)
      
      all_results[[length(all_results) + 1]] <- df_long
    }
  }
  
  do.call(rbind, all_results)
}
df <- process_ecotyper_nested("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/immune_cell_infiltration")

table(df$Sample)
table(df$Sample, df$CellType) # row counts

ggplot(df, aes(x = CellType, y = Abundance)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#BOXPLOT OF CELL STATES
cell_types_of_interest <- c(
  "Monocytes.and.Macrophages",
  "NK.cells",
  "Dendritic.cells",
  "CD4.T.cells",
  "CD8.T.cells"
)

plot_data <- df[df$CellType %in% cell_types_of_interest, ]
plot_data <- plot_data[!is.na(plot_data$Abundance) & !is.na(plot_data$Cell.State), ]
plot_data$Cell.State <- as.factor(plot_data$Cell.State)
# combined OCCC
cell_states <- sort(unique(plot_data$Cell.State))
cell_state_colors <- setNames(brewer.pal(n = length(cell_states), name = "Paired"), cell_states)

# violin plot
ggplot(plot_data, aes(x = Cell.State, y = Abundance, fill = Cell.State)) +
  geom_violin(trim = FALSE) + geom_boxplot(width = 0.1, outlier.size = 0.3) +
  facet_wrap(~ CellType, scales = "free_x", ncol = 2) + labs(x = "Cell State", y = "Abundance") +
  scale_fill_manual(values = cell_state_colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x= element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(1, "lines"))



# FOR EACH 5 OCCC datasets
# colour by cell states
cell_states <- sort(unique(plot_data$Cell.State))
cell_state_colors <- setNames(brewer.pal(n = length(cell_states), name = "Paired"), cell_states)
# map the names 
sample_labels <- c(
  "ecotyper_output_CCOC" = "Bolton et al. (2022)",
  "ecotyper_output_E-MTAB_RNAseq_mRNA" = "Expression Atlas (2026)",
  "ecotyper_output_GSE160692" = "GSE160692",
  "ecotyper_output_GSE189553" = "GSE189553",
  "ecotyper_output_SD1"= "Nagasawa et al. (2019)")

plot_data$Study <- sample_labels[plot_data$Sample]
plot_data$Study <- factor(plot_data$Study, levels = names(ann_colors$Study))

ggplot(plot_data, aes(x = Cell.State, y = Abundance, fill = Cell.State)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4) +
  facet_grid(Study ~ CellType, scales = "free_x") + labs(x= "Cell State", y= "Abundance"
  ) + scale_fill_manual(values = cell_state_colors) + theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text.x = element_text(face = "bold", size = 9),
    strip.text.y = element_text(face = "bold", size = 8, angle = 0),
    panel.spacing = unit(0.5, "lines"))


###################

# ECOTYPE



process_ecotypes <- function(root_path) {
  
  samples <- list.files(root_path, full.names = TRUE)
  all_results <- list()
  
  for (sample_path in samples) {
    
    ecotype_path <- file.path(sample_path, "Carcinoma_Ecotypes")
    if (!dir.exists(ecotype_path)) next
    
    abd_file <- list.files(ecotype_path, pattern = "Abundance", full.names = TRUE)
    assn_file <- list.files(ecotype_path, pattern = "Assignment", full.names = TRUE)
    
    if (length(abd_file) == 0 || length(assn_file) == 0) next
    
    abd <- tryCatch(read.delim(abd_file[1], check.names = FALSE),
                    error = function(e) NULL)
    
    assn <- tryCatch(read.delim(assn_file[1], check.names = FALSE),
                     error = function(e) NULL)
    
    if (is.null(abd) || is.null(assn)) next
    
    colnames(abd)[1] <- "ID"
    colnames(assn)[1] <- "ID"
    
    colnames(assn) <- c("ID", "Ecotype")

    abd_long <- melt(
      abd,
      id.vars = "ID",
      variable.name = "Ecotype_Component",  # CE1–CE10
      value.name = "Abundance"
    )

    df_long <- merge(abd_long, assn, by = "ID", all.x = TRUE)
    
    df_long$Sample <- basename(sample_path)
    
    all_results[[length(all_results) + 1]] <- df_long
  }
  
  do.call(rbind, all_results)
}

ecotype_df <- process_ecotypes("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/immune_cell_infiltration")





########
# ecotypes enriched in some?
assn_all <- unique(ecotype_df[, c("ID", "Sample", "Ecotype")]) # ecotype per sample
table(assn_all$Ecotype)
table(assn_all$Ecotype, assn_all$Sample)
prop.table(table(assn_all$Ecotype, assn_all$Sample), margin=2)


sample_labels <- c(
  "ecotyper_output_CCOC" = "Bolton et al. (2022)",
  "ecotyper_output_E-MTAB_RNAseq_mRNA" = "Expression Atlas (2026)",
  "ecotyper_output_GSE160692" = "GSE160692",
  "ecotyper_output_GSE189553" = "GSE189553",
  "ecotyper_output_SD1"= "Nagasawa et al. (2019)")

assn_all$Sample_label <- sample_labels[assn_all$Sample]

plot_df <- assn_all %>% count(Sample_label, Ecotype) %>%
  group_by(Sample_label) %>% mutate(prop = n / sum(n))

ggplot(plot_df, aes(x = Sample_label, y = prop, fill = Ecotype)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.3) +
  geom_text(aes(label = scales::percent(prop)),
            position = position_stack(vjust = 0.5),
            size = 3) +
  scale_fill_brewer(palette = "Set3") +
  theme_bw() +
  labs(
    x = "Dataset",
    y = "Proportion of samples",
    fill = "Ecotype"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

# heatmap
sample_labels <- c(
  "ecotyper_output_CCOC" = "Bolton et al. (2022)",
  "ecotyper_output_E-MTAB_RNAseq_mRNA" = "Expression Atlas (2026)",
  "ecotyper_output_GSE160692" = "GSE160692",
  "ecotyper_output_GSE189553" = "GSE189553",
  "ecotyper_output_SD1"= "Nagasawa et al. (2019)")

ecotype_df$Sample_label <- sample_labels[ecotype_df$Sample]

heatmap_df <- ecotype_df %>%
  group_by(Sample_label, Ecotype_Component) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")

heatmap_mat <- heatmap_df %>%
  pivot_wider(names_from = Ecotype_Component,
              values_from = mean_abundance)

rownames_mat <- heatmap_mat$Sample_label
heatmap_mat <- as.data.frame(heatmap_mat[, -1])
rownames(heatmap_mat) <- rownames_mat

pheatmap(as.matrix(heatmap_mat),
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")



###########################
# Pathway connection

ecotype_per_sample <- assn_all %>%
  dplyr::select(ID, Sample, Ecotype) %>%
  mutate(Sample_label = sample_labels[Sample])

# paths on col, sample on row
ssgsea_t <- as.data.frame(t(results_ssgsea))
ssgsea_t$ID <- rownames(ssgsea_t)
# merge
merged_df <- left_join(ecotype_per_sample, ssgsea_t, by = "ID")

hallmark_cols <- colnames(ssgsea_t)[colnames(ssgsea_t) != "ID"]

heatmap_input <- merged_df %>%
  filter(!is.na(Ecotype)) %>% 
  group_by(Ecotype) %>%
  summarise(across(all_of(hallmark_cols), mean, na.rm = TRUE)) %>%
  column_to_rownames("Ecotype")


pheatmap(t(heatmap_input),
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         fontsize_row = 7)



table(is.na(merged_df$Ecotype))
# which samples na
merged_df %>% filter(is.na(Ecotype)) %>% distinct(ID, Sample_label)

head(rownames(ssgsea_t))
head(ecotype_per_sample$ID)
sum(ecotype_per_sample$ID %in% rownames(ssgsea_t))
sum(rownames(ssgsea_t) %in% ecotype_per_sample$ID)
intersect(ecotype_per_sample$ID, rownames(ssgsea_t)) %>% length()

###########################
# abundance weighted proportions of cell types for ecotypes 

merged_ct <- df %>%
  dplyr::select(ID, CellType, Abundance) %>%
  left_join(assn_all[, c("ID", "Ecotype")], by = "ID")

plot_df <- merged_ct %>%
  filter(!is.na(Ecotype)) %>%
  group_by(Ecotype, CellType) %>%
  summarise(total_abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Ecotype) %>%
  mutate(prop = total_abundance / sum(total_abundance))


ggplot(plot_df, aes(x = Ecotype, y = prop, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_classic() +
  labs(x = "Ecotype", y = "Proportion of cells", fill = "Cell Type") +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(plot_df, aes(x = Ecotype, y = prop, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.7, colour = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Paired") +
  theme_classic() +
  labs(
    x = "Ecotype",
    y = "Proportion of cells",
    fill = "Cell Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

###################
# exploration

pathway_of_interest <- "HALLMARK_KRAS_SIGNALING_DN"


ggplot(merged_df, aes(x = Ecotype, y = .data[[pathway_of_interest]], fill = Ecotype)) +
  geom_boxplot() +
  facet_wrap(~ Sample_label) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

