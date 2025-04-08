#' Compare workflow results with Zhao et al
library(data.table)
library(glue)
library(ggplot2)
library(ggforestplot)
library(tidyverse)
library(yaml)

params <- yaml::read_yaml("./parameters/script_params.yml")[['zhao_etal']]

tfile <- file.path(tempdir(),'tmp.png')
outdir <- params$outdir

# standardized mask naming convention
zhao_mask_info <- c(
  "HC_PTV" = "pLoF (AAF < 0.1%)",
  "MISS_REVEL0_7" = "REVEL > 0.7 (AAF < 0.1%)",
  "MISS_REVEL0_5" = "REVEL > 0.5 (AAF < 0.1%)"
)
pipeline_mask_info <- c(
  "zhao_plof.0.001" = "pLoF (AAF < 0.1%)",
  "zhao_revel_07.0.001" = "REVEL > 0.7 (AAF < 0.1%)",
  "zhao_revel_05.0.001" = "REVEL > 0.5 (AAF < 0.1%)"
)
zhao_results <- read.table(
  params$zhao_res,sep="\t",header=T,quote=""
)
zhao_results$matched_masks <- zhao_mask_info[zhao_results$MASK]

# obtain the corresponding entry as ExWAS ran with multiple AAF filter
# # e.g., singletons
pipeline_res_dir <- params$pipeline_res
shared_data <- data.frame()
for (each_r in 1:nrow(zhao_results)){
  res <-  zhao_results[each_r,]
  gene <- res$ENSG
  if (res$matched_masks == "pLoF (AAF < 0.1%)"){
    pipeline_file <- read.table(
      file.path(pipeline_res_dir,"REGENIE_OUTPUTS","Regenie_S2","zhao_plof",glue("8_regenie_S2_OUT_wes_qc_chr{res$chrom}_BMI.regenie.gz")),
      sep="\t",header=T,quote=""
    )
    results <- pipeline_file %>% filter(
      Name == glue("{gene}.zhao_plof.0.001")
    )
  }else if (res$matched_masks == "REVEL > 0.7 (AAF < 0.1%)"){
    pipeline_file <- read.table(
      file.path(pipeline_res_dir,"REGENIE_OUTPUTS","Regenie_S2","zhao_revel_07",glue("8_regenie_S2_OUT_wes_qc_chr{res$chrom}_BMI.regenie.gz")),
      sep="\t",header=T,quote=""
    )
    results <- pipeline_file %>% filter(
      Name == glue("{gene}.zhao_revel_07.0.001")
    )
  }else if (res$matched_masks == 'REVEL > 0.5 (AAF < 0.1%)'){
    pipeline_file <- read.table(
      file.path(pipeline_res_dir,"REGENIE_OUTPUTS","Regenie_S2","zhao_revel_05",glue("8_regenie_S2_OUT_wes_qc_chr{res$chrom}_BMI.regenie.gz")),
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
shared_data$matched_masks = pipeline_mask_info[shared_data$Alt]
shared_data <- shared_data %>% dplyr::select(
  ENSG,symbol,matched_masks,Effect,Pval,pipe_se,AAF
) %>% dplyr::rename(
  SE = pipe_se
) %>% dplyr::mutate(
  source = 'Pipeline'
)

# standardized column names
zhao_results <- zhao_results %>% dplyr::select(
  ENSG,SYMBOL,matched_masks,BETA..kg.m2.,SE,P_BOLT_LMM_INF,A1FREQ
) %>% dplyr::rename(
  Effect = BETA..kg.m2.,
  Pval = P_BOLT_LMM_INF,
  symbol = SYMBOL,
  AAF = A1FREQ
) %>% dplyr::mutate(
  source = 'Zhao et al 2024'
)
p_data <- rbind(
  shared_data,
  zhao_results
)
beta_forest_plots <- ggforestplot::forestplot(
  df = p_data,
  name = matched_masks,
  estimate = Effect,
  se = SE,
  pvalue = Pval,
  psignif = 0.05,
  colour = source
) + facet_wrap(~symbol) + theme_classic() +
theme(text = element_text(size=20)) +
labs(
  x = "ExWAS gene-burden test effect size"
)
ggsave(
  tfile,
  beta_forest_plots,
  device='png',width=10,height=7,
)

ggsave(
  file.path(outdir,"zhao_etal_replication.pdf"),
  beta_forest_plots,
  device='pdf',width=10,height=8,
)


# comparing p-values
pval_pval_plots <- inner_join(
  shared_data,
  zhao_results,
  by = c("ENSG","matched_masks")
)
pval_pval_plots <- pval_pval_plots %>% 
  dplyr::mutate(
    nlog10p_pipeline = -1 * log10(Pval.x),
    nlog10p_zhao = -1 * log10(Pval.y)
  )
pearson_r2 <- cor(
  pval_pval_plots$nlog10p_pipeline,
  pval_pval_plots$nlog10p_zhao
)
p <- ggplot(pval_pval_plots) +
  geom_point(
    aes(x=nlog10p_pipeline,y=nlog10p_zhao)
  ) + theme_classic() + geom_abline(slope=1,intercept=0,linetype=2) +
  xlim(0,25) +
  ylim(0,25)+
  theme(text = element_text(size=20)) +
  labs(
    x = bquote("-log"[10]~"(pipeline p-value)"),
    y = bquote("-log"[10]~"(Zhao et al 2024 p-value)"),
    title = glue("Pearson correlation {round(pearson_r2,2)}")
  )
ggsave(
  tfile,p,
  width=8,height=6
)


ggsave(
  file.path(outdir,'zhao_etal_replication_log10p.pdf'),
  p,
  device='pdf',
  width=8,height=6
)
