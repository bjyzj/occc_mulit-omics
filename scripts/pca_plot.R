

library(tidyverse)
library(ggplot2)


sample_info <- readRDS("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/sample_info_OCandRC.rds")
expr_mat <- readRDS("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/expression_matrix_occc_renal.rds")

sample_info_ovary_only <- readRDS("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/sample_info_occc_ovary.rds")
expr_mat_ovary_only <- readRDS("/Users/beyzaerkal/Desktop/occc_multi-omics/processed/expression_matrix_occc_ovary.rds")

pca <- prcomp(t(expr_mat), scale. = TRUE)

plot(pca, type = "l") # scree plot
summary(pca)
summary(pca)$importance[2, 1:2] # variance proportion pc1 and pc2 = (pca$sdev^2) / sum(pca$sdev^2)
# top genes show transcriptional differences
head(pca$rotation[,1])
head(pca$rotation[,2])

pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], sample_info)

ggplot(pca_df, aes(PC1, PC2, color = Tissue, shape = SourceType)) +
  geom_point(size = 3) +
  theme_minimal()

# Assuming multivariate normal distribution for the ellipses. 
# change the label name
mapping <- c("GTEx_Ovary" = "Normal Ovary", "GTEx_Renal_Cortex" = "Normal Renal Cortex")
sample_info <- sample_info %>% mutate(across(c(Study, Group), ~ recode(.x, !!!mapping)))
# 4 groups 2 normal 2 cancer
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], sample_info)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  stat_ellipse(aes(group = Group), type = "norm", level = 0.95, alpha = 0.1, geom = "polygon") + 
  geom_point(aes(shape = Study), size = 2, alpha = 0.7) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 15, 16, 17)) + 
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(shape = "Study Source",
       color = "Sample Group",
       fill = "Sample Group")
