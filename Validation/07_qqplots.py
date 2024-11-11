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



pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/PAST_pipeline_results/validation_10_traits_plof_or_5in5/regeneron/Regenie_S2"
pipeline_result_files = {
  # "ZBMD":"8_regenie_S2_OUT_wes_qc_chr*_ZBMD.regenie.gz",
  # "dBilirubin":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_biliru.regenie.gz",
  # "Calcium":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_Ca.regenie.gz",
  "DBP":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_DBP.regenie.gz"
  # "Glucose":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_glu.regenie.gz",
  # "Height":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_height.regenie.gz",
  # "LDL":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_LDL.regenie.gz",
  # "RBC":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_RBC.regenie.gz",
  # "SBP":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_SBP.regenie.gz",
  # "TG":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_TG.regenie.gz"
}
for trait in pipeline_result_files.keys():
  all_pipeline_results = glob.glob(os.path.join(pipeline_result_path,pipeline_result_files[trait]))
  pipeline_results = pd.DataFrame()
  for each_f in all_pipeline_results:
    each_res = pd.read_csv(
      os.path.join(pipeline_result_path,each_f),
      sep="\t",quoting=csv.QUOTE_NONE,comment='#'
    ).assign(
      LOG10P = lambda df: df['Pval'].apply(lambda val: -1 * np.log10(val))
    )[['Name',"Alt","Model","Chr",'Pos',"Ref","LOG10P","Pval","Effect"]].rename(
      {
        "Name":"Name_pipeline",
        "Alt":"Masks",
        "Chr":"chromosome",
        "Pos":"position_pipeline",
        "Effect":"Beta (Pipeline results)",
        "Model":"Models",
        "Pval":"Pval (Pipeline results)",
        "LOG10P":"Log10P (Pipeline results)"
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
  __makeqq(
    pipeline_results.query(f"Masks == '{m}'"),
    'Pval (Pipeline results)',
    m)
  input()

backman_gwas_path=pd.read_csv("/home/richards/kevin.liang2/scratch/exwas_pipeline/data/backman_exwas_results/GCST90083131_buildGRCh38.tsv.gz",sep="\t",quoting=csv.QUOTE_NONE)
for m in backman_gwas_path.effect_allele.unique():
  __makeqq(
    backman_gwas_path.query(f"effect_allele == '{m}'"),
    'p_value',m
  )
  input()

