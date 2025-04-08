#' Compare results that were no longer significant in workflow but were identified by Backman et al
library(data.table)
library(glue)
library(tidyverse)
library(ggExtra)
library(patchwork)
library(yaml)

params <- yaml::read_yaml("./parameters/script_params.yml")[['compare_nonsig_backman']]


tab_out <- outdir <- params$outdir
tfile <- file.path(tempdir(),'tmp.png')

m3_plot_ofile <- file.path(tab_out,'fig_s5a_backman_inconsistent_compare_m3.pdf')
m1_plot_ofile <- file.path(tab_out,'fig_s5b_backman_inconsistent_compare_m1.pdf')


n_masks <- 2
n_frequencies <- 5
n_burden_tests <- 1
n_phenotypes <- 9
backman_bonferroni <- 0.05/(20000 * n_masks * n_frequencies * n_burden_tests * n_phenotypes)
backman_res <- read.table(
  file.path(outdir,'all_backman_results.tsv.gz'),sep="\t",quote="",
  header=T
) %>% dplyr::filter(Trait != "ZBMD")
pipeline_backman_res <- read.table(
  file.path(outdir,'all_pipeline_backman_res.tsv.gz'),sep="\t",quote="",header=T
)
assertthat::assert_that(
  identical(
    sort(backman_res$Trait %>% unique()),
    sort(pipeline_backman_res$Trait %>% unique())
  )
)
assertthat::assert_that(
  identical(
    sort(backman_res$Models %>% unique()),
    sort(pipeline_backman_res$Models %>% unique())
  )
)
assertthat::assert_that(
  identical(
    sort(backman_res$Masks %>% unique()),
    sort(pipeline_backman_res$Masks %>% unique())
  )
)
backman_sig_plot_data <- data.frame()
n_traits <- sort(pipeline_backman_res$Trait %>% unique())
for (each_trait in n_traits){
  # backman results
  backman_trait_res <- backman_res %>% 
    dplyr::filter(
      Trait == each_trait
    )
  backman_trait_res$Name_backman <- gsub(
    "^.*?\\(|\\).*$","",backman_trait_res$Name_backman
  )
  backman_n_significant <- backman_trait_res %>% dplyr::filter(
    `Pval..Backman.et.al.2021.` < backman_bonferroni
  )
  
  # obtained the corresponding entry in pipeline
  merged_df <- inner_join(
    backman_n_significant,
    pipeline_backman_res,
    by=c("Name_backman" = "Name_pipeline","Masks","Models",'Trait')
  )
  setdiff(backman_n_significant$Name_backman,merged_df$Name_backman)


  backman_sig_plot_data <- rbind(
    backman_sig_plot_data,
    merged_df %>% dplyr::filter(Pval..Pipeline.results. >= backman_bonferroni)
  )
}

m1_p_data <- backman_sig_plot_data %>% filter(grepl("^M1",Masks))
pm1 <- ggplot(m1_p_data) + geom_point(
  aes(
      x = LOG10P..Pipeline.results.,y=LOG10P..Backman.et.al.2021.
    )
  ) +
  geom_vline(xintercept = -1*log10(backman_bonferroni),linetype=2)+
  theme_classic() +
  theme(
    text = element_text(size=20)
  )+
  labs(
    x = bquote("-log"[10]~"(Pipeline p-value)"),
    y = bquote("-log"[10]~"(Backman et al 2021 p-value)")
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

pm3_data <- backman_sig_plot_data %>% filter(grepl("^M3",Masks))
pm3 <- ggplot(pm3_data) + geom_point(
    aes(
      x = LOG10P..Pipeline.results.,
      y=LOG10P..Backman.et.al.2021.
    )
  ) +
  geom_vline(xintercept = -1*log10(backman_bonferroni),linetype=2)+
  theme_classic() +
  theme(
    text = element_text(size=20)
  ) +
  labs(
    x = bquote("-log"[10]~"(Pipeline p-value)"),
    y = bquote("-log"[10]~"(Backman et al 2021 p-value)")
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