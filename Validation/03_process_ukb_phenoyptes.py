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
  p21001_i0: BMI0
  p48_i0: Waist_circumference
  p49_i0: Hip_circumference
"""

import os,sys,shutil,yaml,re,glob,csv,gzip
import pandas as pd
import numpy as np

# inputs
input_file = "/home/richards/kevin.liang2/scratch/exwas_pipeline/data/UKB_phenotypes/UKB_phenotype/UKB_continuous_trait_Nov232023_participant.csv"
# outputs
output_dir = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/processed_UKB_phenotypes"



raw_phenotypes = pd.read_csv(input_file,sep=",",header=0,quoting=csv.QUOTE_NONE)

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

rel_phenotypes = raw_phenotypes[["FID","IID",'SBP','DBP','Standing_height','LDL',"TG","Calcium","Dbilirubin","Glucose","RBC","BMI"]]

with gzip.open(os.path.join(output_dir,'UKB_phenotypes_renamed_columns.tsv.gz'),'wt') as ptr:
  rel_phenotypes.to_csv(
    ptr,quoting=csv.QUOTE_NONE,header=True,index=False,sep="\t",na_rep="NA"
  )