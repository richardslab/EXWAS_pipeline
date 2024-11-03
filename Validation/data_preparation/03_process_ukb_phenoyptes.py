"""
  Process the input phenotypes.

  obtained from Yiheng's Alphamissense paper. The phenoyptes are:

  p4080_i0_a0: SBP
  p4079_i0_a0: DBP
  p50_i0: Standing_height
  p30780_i0: LDL
  p30870_i0: TG
  p30680_i0: Calcium
  p30660_i0: Dbilirubin
  p31: Sex
  p22001: Genetic_sex
  p21022: Age_at_recruitment
  p30740_i0: Glucose
  p30010_i0: RBC
  p21001_i0: BMI
  p48_i0: Waist_circumference
  p49_i0: Hip_circumference

  Processing: 
    - IRNT transform the data.
    - Among EUR
"""

import os,sys,shutil,yaml,re,glob,csv,gzip
import pandas as pd
import numpy as np
from scipy import stats

# inputs
input_file = "/home/richards/kevin.liang2/scratch/exwas_pipeline/data/UKB_phenotypes/UKB_phenotype/UKB_continuous_trait_Nov232023_participant.csv"
eur_individuals="/project/richards/guillaume.butler-laporte/ukb_covid_gwas/anc/ukb.eurFIDIIDPCA.txt"
# outputs
output_dir = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/processed_UKB_phenotypes"

# EUR individuals
eur_ind = pd.read_csv(eur_individuals,sep=" ",header=None,dtype={0:str,1:str}).rename(
  {
    0:"FID",1:"IID"
  },axis=1
)

# Rename the columns
raw_phenotypes = pd.read_csv(input_file,sep=",",header=0,quoting=csv.QUOTE_NONE,dtype={"eid":str})

columns_info = {
  "eid":"FID",
  "p4080_i0_a0": "SBP",
  "p4079_i0_a0": "DBP",
  "p50_i0": "Standing_height",
  "p30780_i0": "LDL",
  "p30870_i0": "TG",
  "p30680_i0": "Calcium",
  "p30660_i0": "Dbilirubin",
  "p31": "Sex",
  "p22001": "Genetic_sex",
  "p21022": "Age_at_recruitment",
  "p30740_i0": "Glucose",
  "p30010_i0": "RBC",
  "p21001_i0": "BMI",
  "p48_i0": "Waist_circumference",
  "p49_i0": "Hip_circumference"
}


raw_phenotypes = raw_phenotypes.rename(
  columns_info,axis=1
)

raw_phenotypes = raw_phenotypes.assign(
  IID=lambda df: df['FID']
)

rel_phenotypes = raw_phenotypes[["FID","IID",'SBP','DBP','Standing_height','LDL',"TG","Calcium","Dbilirubin","Glucose","RBC","BMI"]].query(
  f"FID.isin({eur_ind.FID.to_list()})"
).copy()


# IRNT transform all the columns
phenotype_cols = ['SBP','DBP','Standing_height','LDL',"TG","Calcium","Dbilirubin","Glucose","RBC","BMI"]
def irnt(pheno_series):
  """Performs Inverse Rank Normal Transformation
  
  Based on: https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html#inverse-normal-transformation

  essentially: PPF((rank-0.5)/n).

  Args:
      pheno_series (pandas series): series of value

  Returns:
      numpy array: an array of IRNT values
  """
  val_ranks = pheno_series.rank() # ranks. top rank = largest value
  rank_transformed = (val_ranks-0.5)/(np.sum(~pheno_series.isna()))
  irnt_values = stats.norm.ppf(rank_transformed)
  return irnt_values
irnt_rel_phenotypes = rel_phenotypes.copy()
irnt_rel_phenotypes[phenotype_cols] = irnt_rel_phenotypes[phenotype_cols].apply(
  lambda col: irnt(col),
  axis=0
)



with gzip.open(os.path.join(output_dir,'UKB_phenotypes_renamed_columns_EUR_IRNT.tsv.gz'),'wt') as ptr:
  irnt_rel_phenotypes.to_csv(
    ptr,quoting=csv.QUOTE_NONE,header=True,index=False,sep="\t",na_rep="NA"
  )