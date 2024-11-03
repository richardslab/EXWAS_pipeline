import os,shutil,yaml,re,glob,csv,tempfile
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
tdir = tempfile.TemporaryDirectory()
tfile = os.path.join(tdir.name,'tmp.png')

alphamiss_phenotypes = pd.read_csv("/scratch/richards/yiheng.chen/project14_ExWAS_AlphaMissense/data/UKB_eur_continuous_phenotypes.txt",sep="\t",quoting=csv.QUOTE_NONE)

pipeline_phenotypes = pd.read_csv("/home/richards/kevin.liang2/scratch/exwas_pipeline/results/processed_UKB_phenotypes/UKB_phenotypes_renamed_columns_EUR_IRNT.tsv.gz",sep="\t",quoting=csv.QUOTE_NONE)

shared_data = pd.merge(
  alphamiss_phenotypes[['FID','IRNT_BMI']],
  pipeline_phenotypes[["FID","BMI"]],
  'inner',on=['FID']
)

# the IRNT results seems to be different
fig,ax = plt.subplots()
sns.scatterplot(
  shared_data,
  x='BMI',y='IRNT_BMI',ax=ax
)
fig.savefig(tfile)
plt.close(fig)