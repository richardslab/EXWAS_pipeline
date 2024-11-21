#' phenotypes are normally distributed after transformations
import os,sys,shutil,yaml,re,glob,csv
import matplotlib.pyplot as plt
import seaborn as sns
import tempfile
import pandas as pd
tdir = tempfile.TemporaryDirectory()
tfile = os.path.join(tdir.name,'tmp.png')

eur_individuals="/project/richards/guillaume.butler-laporte/ukb_covid_gwas/anc/ukb.eurFIDIIDPCA.txt"
pheno_file = "/home/richards/kevin.liang2/scratch/exwas_pipeline/data/UKB_phenotypes/UKB_phenotype/UKB_continuous_trait_Nov232023_participant.csv"
eur_ind = pd.read_csv(eur_individuals,sep=" ",header=None,dtype={0:str,1:str}).rename(
  {
    0:"FID",1:"IID"
  },axis=1
)

# Rename the columns
raw_phenotypes = pd.read_csv(pheno_file,sep=",",header=0,quoting=csv.QUOTE_NONE,dtype={"eid":str})

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


phenotypes = ['SBP','DBP','Standing_height','LDL',"TG","Calcium","Dbilirubin","Glucose","RBC","BMI"]
raw_phenotypes = raw_phenotypes.rename(
  columns_info,axis=1
)
raw_phenotypes = raw_phenotypes.assign(
  IID=lambda df: df['FID']
)
rel_phenotypes = raw_phenotypes[["FID","IID"]+phenotypes].query(
  f"FID.isin({eur_ind.FID.to_list()})"
).copy()

b


row_idx = [x%3 for x in range(0,3)]
col_idx = [x%4 for x in range(0,4)]
indices = [(r,c) for r in row_idx for c in col_idx]
fig,ax = plt.subplots(3,4,figsize=(10,10))
for i,t in enumerate(phenotypes):
  sns.histplot(
    rel_phenotypes,
    x = t,ax = ax[indices[i]]
  )

fig.savefig(tfile)

transformed_phenotypes = pd.read_csv("/home/richards/kevin.liang2/scratch/exwas_pipeline/results/processed_UKB_phenotypes/UKB_phenotypes_renamed_columns_EUR_IRNT.tsv.gz",sep="\t",quoting=csv.QUOTE_NONE)

row_idx = [x%3 for x in range(0,3)]
col_idx = [x%4 for x in range(0,4)]
indices = [(r,c) for r in row_idx for c in col_idx]
fig,ax = plt.subplots(3,4,figsize=(10,10))
for i,t in enumerate(phenotypes):
  sns.histplot(
    transformed_phenotypes,
    x = t,ax = ax[indices[i]]
  )

fig.savefig(tfile)
