library(glue)
library(data.table)
library(glue)
library(tidyverse)

train_fam <- read.table("/home/richards/kevin.liang2/scratch/clsa_gwas_required_files/Cleaned_EUR_train.fam",sep=" ",header=F,quote="")
test_fam <- read.table("/home/richards/kevin.liang2/scratch/clsa_gwas_required_files/Cleaned_EUR_test.fam",sep=" ",header=F,quote="")
val_fam <- read.table("/home/richards/kevin.liang2/scratch/clsa_gwas_required_files/Cleaned_EUR_validation.fam",sep=" ",header=F,quote="")

sum(nrow(train_fam),nrow(test_fam),nrow(val_fam))
nrow(train_fam)/sum(nrow(train_fam),nrow(test_fam),nrow(val_fam))
nrow(test_fam)/sum(nrow(train_fam),nrow(test_fam),nrow(val_fam))
nrow(val_fam)/sum(nrow(train_fam),nrow(test_fam),nrow(val_fam))


clsa_anno <- read.table("/scratch/richards/yiheng.chen/CLSA_metabolomics_data/Aug2024_metabolomics_v2_data_check/CLSA_meta_data_with_failed_batch_annot_CLSA_v1_data.csv",sep=",",header=T,quote="")
unique(clsa_anno$X.fialed_batch.)
problem_batches <- clsa_anno %>% filter(X.fialed_batch. == '"yes"')

clsa_met_info <- read.table("/scratch/richards/kevin.liang2/clsa_gwas_required_files/clsa_plink_sample_id_map.csv",sep=",",header=T,quote="")
clsa_id_map <- setNames(
  clsa_met_info$ADM_GWAS_COM,nm=clsa_met_info$PARENT_SAMPLE_NAME
)

clsa_met <- read.table("/scratch/richards/kevin.liang2/clsa_gwas_required_files/clsa_metabolite_measurements.csv",sep=",",header=T,quote="")

clsa_met$id <- clsa_id_map[clsa_met$PARENT_SAMPLE_NAME]
sum(is.na(clsa_met$id))
clsa_met %>% filter(is.na(id)) %>% pull(PARENT_SAMPLE_NAME)

problem_batches$X.. <- unlist(
  lapply(
    problem_batches$X..,
    FUN = function(x){as.numeric(gsub('"',"",x))}
  )
)



length(intersect(train_fam$V1,problem_batches$X.ADM_GWAS3_COM.))/nrow(train_fam)
intersect(test_fam$id,problem_batches$X.ADM_GWAS3_COM.)
intersect(val_fam$id,problem_batches$X.ADM_GWAS3_COM.)


meta_info <- read.table("/project/richards/restricted/clsa/2006016/2006016_McGillU_BRichards_Baseline/2006016_McGillU_BRichards_CoP6_Baseline.csv",sep=",",header=T)

