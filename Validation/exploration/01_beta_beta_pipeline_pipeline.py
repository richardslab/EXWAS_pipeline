"""
  Within same pipeline, replicatable
"""
#%%

import os,shutil,yaml,pickle,re,glob,csv,gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tempfile

tempdir = tempfile.TemporaryDirectory()
tfile = os.path.join(tempdir.name,'temp.png')



# own regenie results:

pipeline_result_files = {
  "BMI":"8_regenie_S2_OUT_wes_qc_chr*_BMI.regenie.gz"
}

# Compare results
plot_data = pd.DataFrame()
for trait in pipeline_result_files.keys():
  all_pipeline_results = glob.glob(os.path.join(pipeline_result_path,pipeline_result_files[trait]))
  pipeline_results_A = pd.DataFrame()
  for each_f in all_pipeline_results:
    each_res = pd.read_csv(
      os.path.join("/home/richards/kevin.liang2/scratch/exwas_pipeline/results/pipeline_results/regeneron/Regenie_S2",each_f),
      sep="\t",quoting=csv.QUOTE_NONE,comment='#'
    )[['Name',"Alt","Model","Chr",'Pos',"Effect","Pval"]].rename(
      {
        "Name":"Name_pipeline",
        "Alt":"Masks",
        "Chr":"chromosome",
        "Pos":"position_pipeline",
        "Model":"Models",
        "Effect":"Beta (Pipeline results)",
        "Pval":"Pval (Pipeline results)"
      },axis=1
    )
    each_res["Name_pipeline"] = each_res["Name_pipeline"].apply(
      lambda val: val.split(".")[0].strip()
    )
    pipeline_results_A = pd.concat(
      [
        pipeline_results_A,
        each_res
      ],axis=0
    ).reset_index(drop=True)
  # b
  all_pipeline_results = glob.glob(os.path.join(pipeline_result_path,pipeline_result_files[trait]))
  pipeline_results_B = pd.DataFrame()
  for each_f in all_pipeline_results:
    each_res = pd.read_csv(
      os.path.join("/home/richards/kevin.liang2/scratch/exwas_pipeline/results/PAST_pipeline_results/Regenie_S2",each_f),
      sep="\t",quoting=csv.QUOTE_NONE,comment='#'
    )[['Name',"Alt","Model","Chr",'Pos',"Effect","Pval"]].rename(
      {
        "Name":"Name_pipeline",
        "Alt":"Masks",
        "Chr":"chromosome",
        "Pos":"position_pipeline",
        "Model":"Models",
        "Effect":"Beta (Pipeline results)",
        "Pval":"Pval (Pipeline results)"
      },axis=1
    )
    each_res["Name_pipeline"] = each_res["Name_pipeline"].apply(
      lambda val: val.split(".")[0].strip()
    )
    pipeline_results_B = pd.concat(
      [
        pipeline_results_B,
        each_res
      ],axis=0
    ).reset_index(drop=True)
  shared_data = pd.merge(
    pipeline_results_A,
    pipeline_results_B,
    'inner',
    left_on = ['Masks','chromosome','Models','Name_pipeline'],
    right_on=['Masks','chromosome','Models','Name_pipeline']
  )
  shared_data['Trait'] = trait
  plot_data = pd.concat(
    [
      plot_data,
      shared_data
    ],axis=0
  ).reset_index(drop=True)
  

row_idx = [x%3 for x in list(range(0,3))]
col_idx = [x%4 for x in list(range(0,4))]
indicies = [(r,c) for r in row_idx for c in col_idx]
fig,ax = plt.subplots(3,4,figsize=(15,15))
for i,t in enumerate(backman_gwas_files.keys()):
  sns.scatterplot(
    plot_data.query(f"Trait == '{t}'"),
    x='Beta (Pipeline results)_x',
    y='Beta (Pipeline results)_y',ax=ax[indicies[i]]
  )
  ax[indicies[i]].axline(
    (0,0),slope = 1,linestyle='--'
  )
  ax[indicies[i]].update(
    {
      "xlabel":"",
      "ylabel":"",
      "title":t
    }
  )

fig.supxlabel("Backman et al 2021")
fig.supylabel("Pipeline results")
fig.savefig(tfile)
