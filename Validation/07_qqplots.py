"""
  Compare the effect sizes my runs and Backman et al 2021
"""
#%%

import os,shutil,yaml,pickle,re,glob,csv,gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import tempfile

tempdir = tempfile.TemporaryDirectory()
tfile = os.path.join(tempdir.name,'temp.png')



pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/pipeline_results/regeneron/Regenie_S2"
pipeline_result_files = {
  "BMI":"8_regenie_S2_OUT_wes_qc_chr*_BMI.regenie.gz"
}
for trait in pipeline_result_files.keys():
  all_pipeline_results = glob.glob(os.path.join(pipeline_result_path,pipeline_result_files[trait]))
  pipeline_results = pd.DataFrame()
  for each_f in all_pipeline_results:
    each_res = pd.read_csv(
      os.path.join(pipeline_result_path,each_f),
      sep=" ",quoting=csv.QUOTE_NONE,comment='#'
    ).query("TEST == 'ADD'").assign(
      Models = lambda df: df['TEST'].apply(lambda val: f"{val}-WGR-LR"),
      Pval = lambda df: df['LOG10P'].apply(lambda val: 10**(-1 * val))
    )[['ID',"ALLELE1","Models","CHROM",'GENPOS',"ALLELE0","LOG10P","Pval","BETA"]].rename(
      {
        "ID":"Name_pipeline",
        "ALLELE1":"Masks",
        "CHROM":"chromosome",
        "GENPOS":"position_pipeline",
        "BETA":"Beta (Pipeline results)",
        "Pval":"Pval (Pipeline results)"
      },axis=1
    )
    each_res["Name_pipeline"] = each_res["Name_pipeline"].apply(
      lambda val: val.split(".")[0].strip()
    )
    each_res['Trait'] = trait
    pipeline_results = pd.concat(
      [
        pipeline_results,
        each_res
      ],axis=0
    ).reset_index(drop=True)

def __makeqq(pvalues,pval_col,m):
  # uniform distribution PPF = proportion
  pvalues = pvalues.assign(
    pval_ranks = lambda df: df[pval_col].rank()/len(df)
  )
  pvalues = pvalues.assign(
    expected_p = lambda df: df['pval_ranks'].apply(lambda val: -1 * np.log10(stats.uniform.ppf(val))),
    observed_p = lambda df: -1*np.log10(df[pval_col])
  )
  fig,ax = plt.subplots()
  sns.scatterplot(
    pvalues,
    x="expected_p",y='observed_p',ax=ax
  )
  ax.axline((0,0),slope=1)
  ax.update({'title':m})
  fig.savefig(tfile)
  plt.close(fig)
  return

# make qq plots
for m in pipeline_results.Masks.unique():
  pvalues = pipeline_results.query(f"Masks == '{m}'")['Pval (Pipeline results)'].to_frame()
  __makeqq(pvalues,'Pval (Pipeline results)',m)
  input()

backman_gwas_path=pd.read_csv("/home/richards/kevin.liang2/scratch/exwas_pipeline/data/backman_exwas_results/GCST90082670_buildGRCh38.tsv.gz",sep="\t",quoting=csv.QUOTE_NONE)
for m in backman_gwas_path.effect_allele.unique():
  __makeqq(
    backman_gwas_path.query(f"effect_allele == '{m}'"),
    'p_value',m
  )
  input()

