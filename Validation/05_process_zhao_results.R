#' comparison of significant results from zhao et al
library(data.table)
library(glue)
library(ggplot2)
library(ggforestplot)
library(tidyverse)

tfile <- file.path(tempdir(),'tmp.png')
mask_info <- c(
  "HC_PTV" = "zhao_plof.0.001",
  "MISS_REVEL0_7" = "zhao_revel_07.0.001",
  "MISS_REVEL0_5" = "zhao_revel_05.0.001"
)
zhao_results <- read.table(
  "/home/richards/kevin.liang2/scratch/exwas_pipeline/data/zhao_data/S1_zhao_results.tsv",sep="\t",header=T,quote=""
)
zhao_results$matched_masks <- mask_info[zhao_results$MASK]

pipeline_res_dir <- "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/zhao_etal_BSN_BMI"
shared_data <- data.frame()
for (each_r in 1:nrow(zhao_results)){
  res <-  zhao_results[each_r,]
  gene <- res$ENSG
  if (res$matched_masks == "zhao_plof.0.001"){
    pipeline_file <- read.table(
      file.path(pipeline_res_dir,"zhao_plof","Regenie_S2",glue("8_regenie_S2_OUT_wes_qc_chr{res$chrom}_BMI.regenie.gz")),
      sep="\t",header=T,quote=""
    )
    results <- pipeline_file %>% filter(
      Name == glue("{gene}.zhao_plof.0.001")
    )
  }else if (res$matched_masks == "zhao_revel_07.0.001"){
    pipeline_file <- read.table(
      file.path(pipeline_res_dir,"zhao_revel_07","Regenie_S2",glue("8_regenie_S2_OUT_wes_qc_chr{res$chrom}_BMI.regenie.gz")),
      sep="\t",header=T,quote=""
    )
    results <- pipeline_file %>% filter(
      Name == glue("{gene}.zhao_revel_07.0.001")
    )
  }else if (res$matched_masks == 'zhao_revel_05.0.001'){
    pipeline_file <- read.table(
      file.path(pipeline_res_dir,"zhao_revel_05","Regenie_S2",glue("8_regenie_S2_OUT_wes_qc_chr{res$chrom}_BMI.regenie.gz")),
      sep="\t",header=T,quote=""
    )
    results <- pipeline_file %>% filter(
      Name == glue("{gene}.zhao_revel_05.0.001")
    )
  }
  shared_data <- rbind(
    shared_data,
    results
  )
}

gene_names <- setNames(zhao_results$SYMBOL,nm=zhao_results$ENSG)

shared_data$ENSG <- gsub("\\..*$","",shared_data$Name)
shared_data$symbol <-  gene_names[shared_data$ENSG]
shared_data$pipe_se = unlist(lapply(
  strsplit(shared_data$Info,";"),
  FUN = function(x){as.numeric(gsub("REGENIE_SE=","",x[[1]]))}
))
shared_data <- shared_data %>% dplyr::select(
  ENSG,symbol,Alt,Effect,Pval,pipe_se,AAF
) %>% dplyr::rename(
  SE = pipe_se
) %>% dplyr::mutate(
  source = 'Pipeline'
)

zhao_results <- zhao_results %>% dplyr::select(
  ENSG,SYMBOL,matched_masks,BETA..kg.m2.,SE,P_BOLT_LMM_INF,A1FREQ
) %>% dplyr::rename(
  Effect = BETA..kg.m2.,
  Pval = P_BOLT_LMM_INF,
  Alt = matched_masks,
  symbol = SYMBOL,
  AAF = A1FREQ
) %>% dplyr::mutate(
  source = 'Zhao et al 2024'
)
p_data <- rbind(
  shared_data,
  zhao_results
)
p <- ggforestplot::forestplot(
  df = p_data,
  name = Alt,
  estimate = Effect,
  se = SE,
  pvalue = Pval,
  psignif = 0.05,
  colour = source
) + facet_wrap(~symbol) + theme_classic() +
theme(text = element_text(size=20))
ggsave(
  tfile,
  p,
  device='png',width=10,height=8,
)

ggsave(
  "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/figures/zhao_etal_replication.png",
  p,
  device='png',width=10,height=8,
)


# comparing p-values
pval_pval_plots <- inner_join(
  shared_data,
  zhao_results,
  by = c("ENSG","Alt")
)
pval_pval_plots <- pval_pval_plots %>% 
  dplyr::mutate(
    nlog10p_pipeline = -1 * log10(Pval.x),
    nlog10p_zhao = -1 * log10(Pval.y)
  )
p <- ggplot(pval_pval_plots) +
  geom_point(
    aes(x=nlog10p_pipeline,y=nlog10p_zhao)
  ) + theme_classic() + geom_abline(slope=1,intercept=0,linetype=2) +
  theme(text = element_text(size=20)) +
  labs(
    x = bquote("-log"[10]~"(pipeline p-value)"),
    y = bquote("-log"[10]~"(Zhao et al 2024 p-value)")
  )
ggsave(
  tfile,p
)
# ggsave(
#   "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/figures/zhao_etal_replication_log10p.png",
#   p,
#   device='png'
# )

# # comparing aaf
# p <- ggplot(pval_pval_plots) +
#   geom_point(
#     aes(x=AAF.x,y=AAF.y)
#   ) + theme_classic() + geom_abline(slope=1,intercept=0,linetype=2) +
#   theme(text = element_text(size=20)) +
#   labs(
#     x = bquote("Pipeline AAF"),
#     y = bquote("Zhao et al AAF")
#   )
# ggsave(
#   tfile,p
# )