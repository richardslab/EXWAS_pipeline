#' Venn diagram of various annotation results
library(data.table)
library(glue)
library(ggplot2)
library(ggVennDiagram)
library(tidyverse)

tfile <- file.path(tempdir(),'tmp.png')
gene_venn <- "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/figures/gene_venn_diagram.png"
pipeline_res <- read.table(
  "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/tables/pipeline_results.tsv.gz",sep="\t",header=T,quote=""
)
alpha_res <- read.table(
  "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/tables/alphamissense_results.tsv.gz",sep="\t",header=T,quote=""
)
regeneron_res <- read.table(
  "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/tables/backman_results.tsv.gz",sep="\t",header=T,quote=""
)

gene_venn_data <- list(
  "Pipeline" = unique(pipeline_res$Name_pipeline),
  "Backman et al\n2021" = unique(regeneron_res$Name_backman),
  "Alphamissense" = unique(alpha_res$Name_Alpha)
)

p <- ggVennDiagram(
  gene_venn_data, 
  color = "black", 
  lwd = 0.8, 
  lty = 1,
  label_alpha = 0.5
) + 
scale_fill_gradient(
  low = "#F4FAFE", 
  high = "#4981BF"
) +
theme_classic()+
theme(
  line = element_blank(),
  text = element_blank(),
  title = element_blank(),
  legend.position = "none"
)
ggsave(
  tfile,
  p,
  device='png',width=7,height=7,units='in'
)
# ggsave(
#   gene_venn,
#   p,
#   device='png',width=7,height=7,units='in'
# )

