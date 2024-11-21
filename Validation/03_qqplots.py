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



pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/plof_or_5in5/regeneron/Regenie_S2"
pipeline_result_files = {
  "ZBMD":"8_regenie_S2_OUT_wes_qc_chr*_ZBMD.regenie.gz",
  "dBilirubin":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_biliru.regenie.gz",
  "Calcium":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_Ca.regenie.gz",
  "DBP":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_DBP.regenie.gz",
  "Glucose":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_glu.regenie.gz",
  "Height":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_height.regenie.gz",
  "LDL":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_LDL.regenie.gz",
  "RBC":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_RBC.regenie.gz",
  "SBP":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_SBP.regenie.gz",
  "TG":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_TG.regenie.gz",
  "WHR":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_WHR.regenie.gz"
}

pipeline_results = pd.DataFrame()
for trait in pipeline_result_files.keys():
  all_pipeline_results = glob.glob(os.path.join(pipeline_result_path,pipeline_result_files[trait]))
  pipeline_results = pd.DataFrame()
  for each_f in all_pipeline_results:
    each_res = pd.read_csv(
      os.path.join(pipeline_result_path,each_f),
      sep="\t",quoting=csv.QUOTE_NONE,comment='#'
    ).query("~Effect.isna()").assign(
      LOG10P = lambda df: df['Pval'].apply(lambda val: -1 * np.log10(val)),
      SE = lambda df: df['Info'].apply(
        lambda val: _parse_info(val)
      )
    )[['Name',"Alt","Model","Chr",'Pos',"Ref","Pval","Effect",'LOG10P','SE']].rename(
      {
        "Name":"Name_pipeline",
        "Alt":"Masks",
        "Chr":"chromosome",
        "Pos":"position_pipeline",
        "Model":"Models",
        "Effect":"Beta (Pipeline results)",
        "Pval":"Pval (Pipeline results)",
        "LOG10P":"LOG10P (Pipeline results)",
        "SE":"SE (Pipeline results)"
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

def __makeqq(df):
  all_m = df.Masks.unique()
  t = df.Trait.unique()
  row_idx = [i%2 for i in range(0,2)]
  col_idx = [c%3 for c in range(0,3)]
  idx = [tuple([r,c]) for r in row_idx for c in col_idx]
  fig,ax = plt.subplots(2,3,figsize=(30,25))
  for i,each_m in enumerate(all_m):
    m_df = df.query(f"Masks == '{each_m}'")
    pvalues = m_df[['Pval (Pipeline results)','Models']]
    # uniform distribution PPF = proportion
    pvalues = pvalues.assign(
      pval_ranks = lambda df: df['Pval (Pipeline results)'].rank()/len(df)
    )
    pvalues = pvalues.assign(
      expected_p = lambda df: df['pval_ranks'].apply(lambda val: -1 * np.log10(stats.uniform.ppf(val))),
      observed_p = lambda df: -1*np.log10(df['Pval (Pipeline results)'])
    )
    sns.scatterplot(
      pvalues,
      x="expected_p",
      y='observed_p',
      ax=ax[idx[i]],
      hue = "Models"
    )
    ax[idx[i]].axline((0,0),slope=1)
    ax[idx[i]].update({'title':f"{each_m} {t}"})
  fig.savefig(tfile)
  plt.close(fig)
  return

# make qq plots
for t in pipeline_results.Trait.unique():
  __makeqq(
    df = pipeline_results.query(f"Trait == '{t}' & ~`Beta (Pipeline results)`.isna()")
  )
  input()
