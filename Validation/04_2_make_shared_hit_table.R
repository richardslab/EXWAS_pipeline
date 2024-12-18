library(data.table)
library(glue)
library(tidyverse)


outdir <- "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/tables" 

tab_out <- "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/figures"


# make backman table
# bonferroni correction per phenotype
n_masks <- 2
n_frequencies <- 5
backman_bonferroni <- 0.05/(20000*n_masks*n_frequencies)
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
backman_summary_table <- data.frame()
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
  backman_most_sig_gene <- backman_n_significant %>% slice_min(`Pval..Backman.et.al.2021.`) %>% select(Name_backman,Masks)
  b_sig_genes <- c()
  for (i in 1:nrow(backman_most_sig_gene)){
    b_sig_genes <- c(
      b_sig_genes,
      paste0(
        c(backman_most_sig_gene$Name_backman[i],
        backman_most_sig_gene$Masks[i]),collapse="-"
      )
    )
  }
  b_sig_genes <- paste0(b_sig_genes,collapse=";")

  
  # pipeline backman results
  pipeline_trait_res <- pipeline_backman_res %>% 
    dplyr::filter(Trait == each_trait)
  pipeline_n_significant <-  pipeline_trait_res %>% 
    dplyr::filter(
      `Pval..Pipeline.results.` < backman_bonferroni  
    )
  pipeline_most_sig <- pipeline_n_significant  %>% slice_min(
    `Pval..Pipeline.results.`
  ) %>% select(Name_pipeline,Masks)
  p_sig_genes <- c()
  for (i in 1:nrow(pipeline_most_sig)){
    p_sig_genes <- c(
      p_sig_genes,
      paste0(
        c(pipeline_most_sig$Name_pipeline[i],
        pipeline_most_sig$Masks[i]),collapse="-"
      )
    )
  }
  p_sig_genes <- paste0(p_sig_genes,collapse=";")

  # shared hits
  shared_sig_hits <- intersect(
    backman_n_significant %>% pull(Name_backman) %>% unique(),
    unique(pipeline_n_significant$Name_pipeline)
  )
  backman_unique <- setdiff(
    backman_n_significant %>% pull(Name_backman) %>% unique(),
    unique(pipeline_n_significant$Name_pipeline)
  )
  pipeline_unique <- setdiff(
    unique(pipeline_n_significant$Name_pipeline),backman_n_significant %>% pull(Name_backman) %>% unique()
  )
  
  backman_summary_table <- rbind(
    backman_summary_table,
    data.frame(
      Trait = each_trait,
      "N hits study" = length(unique(backman_n_significant$Name_backman)),
      "N hits pipeline" = length(unique(pipeline_n_significant$Name_pipeline)),
      "Shared hits" = length(shared_sig_hits),
      "unique to study" = length(backman_unique),
      "unique to pipeline" = length(pipeline_unique),
      "Strongest signal study" = b_sig_genes,
      "Strongest signal pipeline" = p_sig_genes
    )
  )  
}

# Make Chen table
n_masks <- 2
n_frequencies <- 3
alphamiss_bonferroni <- 0.05/(20000*n_masks*n_frequencies)
alphamiss_res <- read.table(
  file.path(outdir,'all_alphamiss_res.tsv.gz'),sep="\t",quote="",
  header=T
) %>% filter(Masks != "")
pipeline_alpha_res <- read.table(
  file.path(outdir,'all_pipeline_alpha_res.tsv.gz'),sep="\t",quote="",header=T
) %>% filter(
  !grepl("M3_",Masks) &
  !grepl(".1e-05",Masks) &
  !grepl("0.0001",Masks)
)
assertthat::assert_that(
  identical(
    sort(alphamiss_res$Trait %>% unique()),
    sort(pipeline_alpha_res$Trait %>% unique())
  )
)
alphamiss_summary_table <- data.frame()
n_traits <- sort(alphamiss_res$Trait %>% unique())
for (each_trait in n_traits){
  # alphamiss results
  alpha_trait_res <- alphamiss_res %>% 
    dplyr::filter(
      Trait == each_trait
    )
  alpha_trait_res$Name_Alphamissense <- gsub(
    "\\..*$","",alpha_trait_res$Name_Alphamissense
  )
  alphamiss_n_significant <- alpha_trait_res %>% dplyr::filter(
    `Pval..Chen.et.al.2024.` < alphamiss_bonferroni
  )
  alpha_most_sig_gene <- alphamiss_n_significant %>% slice_min(`Pval..Chen.et.al.2024.`) %>% select(Name_Alphamissense,Masks)
  a_sig_genes <- c()
  for (i in 1:nrow(alpha_most_sig_gene)){
    a_sig_genes <- c(
      a_sig_genes,
      paste0(
        c(alpha_most_sig_gene$Name_Alphamissense[i],
        alpha_most_sig_gene$Masks[i]),collapse="-"
      )
    )
  }
  a_sig_genes <- paste0(a_sig_genes,collapse=";")

  
  # pipeline backman results
  pipeline_trait_res <- pipeline_alpha_res %>% 
    dplyr::filter(Trait == each_trait)
  pipeline_n_significant <-  pipeline_trait_res %>% 
    dplyr::filter(
      `Pval..Pipeline.results.` < alphamiss_bonferroni  
    )
  pipeline_most_sig <- pipeline_n_significant  %>% slice_min(
    `Pval..Pipeline.results.`
  ) %>% select(Name_pipeline,Masks)
  p_sig_genes <- c()
  for (i in 1:nrow(pipeline_most_sig)){
    p_sig_genes <- c(
      p_sig_genes,
      paste0(
        c(pipeline_most_sig$Name_pipeline[i],
        pipeline_most_sig$Masks[i]),collapse="-"
      )
    )
  }
  p_sig_genes <- paste0(p_sig_genes,collapse=";")

  # shared sig hits
  shared_sig_hits <- intersect(
    unique(alphamiss_n_significant$Name_Alphamissense),
    unique(pipeline_n_significant$Name_pipeline)
  )
  alphamiss_unique <- setdiff(
    unique(alphamiss_n_significant$Name_Alphamissense),
    unique(pipeline_n_significant$Name_pipeline)
  )
  pipeline_unique <- setdiff(
    unique(pipeline_n_significant$Name_pipeline),
    unique(alphamiss_n_significant$Name_Alphamissense)
  )
  
  alphamiss_summary_table <- rbind(
    alphamiss_summary_table,
    data.frame(
      Trait = each_trait,
      "N hits study" = length(unique(alphamiss_n_significant$Name_Alphamissense)),
      "N hits pipeline" = length(unique(pipeline_n_significant$Name_pipeline)),
      "Shared hits" = length(shared_sig_hits),
      "unique to study" = length(alphamiss_unique),
      "unique to pipeline" = length(pipeline_unique),
      "Strongest signal study" = a_sig_genes,
      "Strongest signal pipeline" = p_sig_genes
    )
  )  
}


# save
write.table(
  backman_summary_table,
  file.path(tab_out,"backman_summary_table.tsv"),
  sep="\t",row.names=F
)
write.table(
  alphamiss_summary_table,
  file.path(tab_out,"alphamiss_summary_table.tsv"),
  sep="\t",row.names=F
)