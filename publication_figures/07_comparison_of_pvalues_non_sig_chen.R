#' Compare results that were no longer significant in workflow but were identified by Chen et al
library(data.table)
library(glue)
library(tidyverse)
library(ggExtra)
library(patchwork)
library(yaml)

params <- yaml::read_yaml("./parameters/script_params.yml")[['compare_nonsig_backman']]


tab_out <- outdir <- params$outdir

tfile <- file.path(tempdir(),'tmp.png')

m3_plot_ofile <- file.path(tab_out,'fig_s6a_chen_inconsistent_compare_m2.pdf')
m1_plot_ofile <- file.path(tab_out,'fig_s6b_chen_inconsistent_compare_m1.pdf')

n_masks <- 2
n_frequencies <- 3
n_burden_tests <- 1
n_phenotypes <- 21
alphamiss_bonferroni <- 0.05/(20000 * n_masks * n_frequencies * n_burden_tests * n_phenotypes)
alphamiss_res <- read.table(
  file.path(outdir,'all_alphamiss_res.tsv.gz'),sep="\t",quote="",
  header=T
) %>% filter(Masks != "" & !is.na(Beta..Chen.et.al.2024.)) 
pipeline_alpha_res <- read.table(
  file.path(outdir,'all_pipeline_alpha_res.tsv.gz'),sep="\t",quote="",header=T
) %>% filter(
  !grepl("M3_",Masks) &
  !grepl(".1e-05",Masks) &
  !grepl("0.0001",Masks) &
  !is.na(Models)
)
assertthat::assert_that(
  identical(
    sort(alphamiss_res$Trait %>% unique()),
    sort(pipeline_alpha_res$Trait %>% unique())
  )
)
assertthat::assert_that(
  identical(
    sort(alphamiss_res$Masks %>% unique()),
    sort(pipeline_alpha_res$Masks %>% unique())
  )
)
assertthat::assert_that(
  identical(
    sort(alphamiss_res$Models %>% unique()),
    sort(pipeline_alpha_res$Models %>% unique())
  )
)
chen_sig_plot_data <- data.frame()
n_traits <- sort(pipeline_alpha_res$Trait %>% unique())
for (each_trait in n_traits){
  # chen results
  alphamiss_trait_res <- alphamiss_res %>% 
    dplyr::filter(
      Trait == each_trait
    )
  alphamiss_trait_res$Name_Alphamissense <- gsub("\\..*$","",alphamiss_trait_res$Name_Alphamissense)
  alpha_n_significant <- alphamiss_trait_res %>% dplyr::filter(
    `Pval..Chen.et.al.2024.` < alphamiss_bonferroni
  )
  
  # obtained the corresponding entry in pipeline
  merged_df <- inner_join(
    alpha_n_significant,
    pipeline_alpha_res,
    by=c("Name_Alphamissense" = "Name_pipeline","Masks","Models",'Trait')
  )
  assertthat::assert_that(nrow(merged_df) == nrow(alpha_n_significant))


  chen_sig_plot_data <- rbind(
    chen_sig_plot_data,
    merged_df %>% dplyr::filter(Pval..Pipeline.results. >= alphamiss_bonferroni)
  )
}

m1_p_data <- chen_sig_plot_data %>% filter(grepl("^M1",Masks))
pm1 <- ggplot(m1_p_data) + geom_point(
  aes(
      x = LOG10P..Pipeline.results.,y=LOG10P..Chen.et.al.2024.
    )
  ) +
  geom_vline(xintercept = -1*log10(alphamiss_bonferroni),linetype=2)+
  theme_classic() +
  theme(
    text = element_text(size=20)
  )+
  labs(
    x = bquote("-log"[10]~"(Pipeline p-value)"),
    y = bquote("-log"[10]~"(Chen et al 2024 p-value)")
  )
pm1 <- ggMarginal(pm1)
ggsave(
  tfile,
  pm1,
  device='png',width=6,height=6
)
ggsave(
  m1_plot_ofile,
  pm1,
  device='pdf',width=6,height=6
)

pm3_data <- chen_sig_plot_data %>% filter(grepl("^M2",Masks))
pm3 <- ggplot(pm3_data) + geom_point(
    aes(
      x = LOG10P..Pipeline.results.,
      y=LOG10P..Chen.et.al.2024.
    )
  ) +
  geom_vline(xintercept = -1*log10(alphamiss_bonferroni),linetype=2)+
  theme_classic() +
  theme(
    text = element_text(size=20)
  ) +
  labs(
    x = bquote("-log"[10]~"(Pipeline p-value)"),
    y = bquote("-log"[10]~"(Chen et al 2024 p-value)")
  )
pm3 <- ggMarginal(pm3)
ggsave(
  tfile,
  pm3,
  device='png',width=6,height=6
)
ggsave(
  m3_plot_ofile,
  pm3,
  device='pdf',width=6,height=6
)