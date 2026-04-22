
library(tidyverse)
library(purrr)
library(ggplot2)
library(readxl)
library(patchwork)

set.seed(123)
# load hallmark results
res_OCvsRC <- read_excel("/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/GSEA_RNA_results_OCvsRC.xlsx", sheet = "Hallmark_MSigDB")
res_OCvsOv <- read_excel("/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/GSEA_OCvsOVARY_RNA_results.xlsx", sheet = "Hallmark_MSigDB")
res_Prot <- read_excel("/Users/beyzaerkal/Desktop/occc_multi-omics/supplementary/GSEA_Prot_OCvsRC_results.xlsx", sheet = "Hallmark_MSigDB")


res_OCvsRC$Comparison <- "RNA-seq: OCCC vs ccRCC"
res_OCvsOv$Comparison <- "RNA-seq: OCCC vs Normal Ovary"
res_Prot$Comparison   <- "Proteomics: OCCC vs ccRCC"



###########
# plot function
###########
plot_hallmark <- function(df, title_text){
  
  df <- as.data.frame(df) %>%
    mutate(
      significance = case_when(
        p.adjust < 0.001 ~ "***",
        p.adjust < 0.01 ~ "**",
        p.adjust < 0.05 ~ "*",
        TRUE ~ ""
      ),
      Description = gsub("HALLMARK_", "", Description),
      Description = gsub("_", " ", Description)
    ) %>%
    arrange(NES)
  
  ggplot(df, aes(x=1, y=reorder(Description, NES), fill=NES)) +
    geom_tile(color="white") +
    geom_text(aes(label=significance), size=4) +
    scale_fill_gradient2(
      low="steelblue",
      mid="white",
      high="firebrick",
      midpoint=0
    ) +
    labs(title=title_text) +
    theme_minimal() +
    theme(
      axis.title=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      panel.grid=element_blank()
    )
}

##########

p1 <- plot_hallmark(res_OCvsRC, "OCCC vs ccRCC\n(RNA-seq)")
p2 <- plot_hallmark(res_OCvsOv, "OCCC vs Normal Ovary\n(RNA-seq)")
p3 <- plot_hallmark(res_Prot,   "OCCC vs ccRCC\n(Proteomics)")

p1 | p2 | p3






