#' examin annotation overlap between Chen et al and workflow
library(data.table)
library(glue)
library(ggplot2)
library(ggVennDiagram)
library(tidyverse)
library(patchwork)
library(yaml)

params <- yaml::read_yaml("./parameters/script_params.yml")[['overlap_var']]

tfile <- file.path(tempdir(),'tmp.png')
var_file <- file.path(params$outdir, 'var_venn_diagram.pdf')

pipeline_plof_or_del5in5 <- read.table(workflow_annot$workflow_annot,sep="\t",header=F,quote="")

alpha_plof_or_del5in5 <- read.table(
  params$chen_annot,
  sep="\t",header=F,quote=""
)

colnames(alpha_plof_or_del5in5) <- paste0("alpha_",c("SNP","GENE","annotation"))
annotation_map <- c(
  "pLoF" = "pLoF",
  "missense.1in5" = "missense.1in5",
  "missense.5in5" = 'deleterious_5_of_5'
)

alpha_plof_or_del5in5['alpha_annotation_matched'] = annotation_map[alpha_plof_or_del5in5$alpha_annotation]

alpha_plof_or_del5in5$id = paste0(
  alpha_plof_or_del5in5$alpha_SNP,"-",
  alpha_plof_or_del5in5$alpha_annotation_matched
)


colnames(pipeline_plof_or_del5in5) = paste0("pipeline_",c("SNP","GENE","annotation"))

pipeline_plof_or_del5in5$id = paste0(
  pipeline_plof_or_del5in5$pipeline_SNP,"-",
  pipeline_plof_or_del5in5$pipeline_annotation
)
  


make_venn <- function(df){
  p <- ggVennDiagram(
    df, 
    color = "black", 
    lwd = 0.8, 
    lty = 1,
    label_alpha = 0.5,
    label_size=5,
    set_size=5
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
  ) + coord_flip()
  return(p)
}

a <- make_venn(
  var_venn_data <- list(
    "pipeline pLoF             " = pipeline_plof_or_del5in5 %>% filter(pipeline_annotation == "pLoF") %>% pull(id),
    "             Chen et al 2024 pLoF" = alpha_plof_or_del5in5 %>% filter(alpha_annotation == "pLoF") %>% pull(id)
  )
)

b <- make_venn(
  var_venn_data <- list(
    "pipeline 5in5             " = pipeline_plof_or_del5in5 %>% filter(pipeline_annotation == "deleterious_5_of_5") %>% pull(id),
    "             Chen et al 2024 5in5" = alpha_plof_or_del5in5 %>% filter(alpha_annotation == "missense.5in5") %>% pull(id)
  )
)


c <- make_venn(
  var_venn_data <- list(
    "pipeline pLoF             " = pipeline_plof_or_del5in5 %>% filter(pipeline_annotation == "pLoF") %>% pull(id),
    "             Chen et al 2024 5in5" = alpha_plof_or_del5in5 %>% filter(alpha_annotation == "missense.5in5") %>% pull(id)
  )
)

d <- make_venn(
  var_venn_data <- list(
    "pipeline 5in5             " = pipeline_plof_or_del5in5 %>% filter(pipeline_annotation == "deleterious_5_of_5") %>% pull(id),
    "             Chen et al 2024 pLoF" = alpha_plof_or_del5in5 %>% filter(alpha_annotation == "pLoF") %>% pull(id)
  )
)

p <- (a|b)/(c|d)

ggsave(
  tfile,
  p,
  device='png',
  height=6,width=8,units='in'
)




ggsave(
  var_file,
  p,
  device='pdf',
  height=6,width=8,units='in'
)

