"""
  Compare the effect sizes my runs and Backman et al 2021
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



#
output_dir = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/PAST_pipeline_results/testruns_chr10_validations"
# downloaded_data_constants
backman_gwas_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/data/backman_exwas_results"
backman_gwas_files = {
  "SBP":"GCST90083132_buildGRCh38.tsv.gz",
  "DBP":"GCST90083131_buildGRCh38.tsv.gz",
  "Standing_height":"GCST90083263_buildGRCh38.tsv.gz",
  "LDL":"GCST90083019_buildGRCh38.tsv.gz",
  "TG":"GCST90083030_buildGRCh38.tsv.gz",
  "Calcium":"GCST90083009_buildGRCh38.tsv.gz",
  "Dbilirubin":"GCST90083007_buildGRCh38.tsv.gz",
  "Glucose":"GCST90083015_buildGRCh38.tsv.gz",
  "RBC":"GCST90082965_buildGRCh38.tsv.gz",
  "BMI":"GCST90082670_buildGRCh38.tsv.gz"
}
# own regenie results:
pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/PAST_pipeline_results/testruns_chr10/regeneron/Regenie_S2"
pipeline_result_files = {
  "SBP":"8_regenie_S2_OUT_wes_qc_chr10_SBP.regenie.gz",
  "DBP":"8_regenie_S2_OUT_wes_qc_chr10_DBP.regenie.gz",
  "Standing_height":"8_regenie_S2_OUT_wes_qc_chr10_Standing_height.regenie.gz",
  "LDL":"8_regenie_S2_OUT_wes_qc_chr10_LDL.regenie.gz",
  "TG":"8_regenie_S2_OUT_wes_qc_chr10_TG.regenie.gz",
  "Calcium":"8_regenie_S2_OUT_wes_qc_chr10_Calcium.regenie.gz",
  "Dbilirubin":"8_regenie_S2_OUT_wes_qc_chr10_Dbilirubin.regenie.gz",
  "Glucose":"8_regenie_S2_OUT_wes_qc_chr10_Glucose.regenie.gz",
  "RBC":"8_regenie_S2_OUT_wes_qc_chr10_RBC.regenie.gz",
  "BMI":"8_regenie_S2_OUT_wes_qc_chr10_BMI.regenie.gz"
}

backman_masks = {
  "M1.singleton" : "M1_LoF.singleton",
  "M1.01" : "M1_LoF.01",
  "M1.001" : "M1_LoF.001",
  "M1.0001" : "M1_LoF.0001",
  'M3.singleton' : "M2_LoF_deleterious.singleton",
  'M3.01' : "M2_LoF_deleterious.01",
  'M3.001' : "M2_LoF_deleterious.001",
  'M3.0001' : "M2_LoF_deleterious.0001"
}

cnames = ["Trait","chromosome",'position_backman','position_pipeline',"Masks","Models","Beta (Backman et al 2021)","Beta (Pipeline results)","Pval (Backman et al 2021)","Pval (Pipeline results)","Name_backman","Name_pipeline"]
# Compare results
if not os.path.isfile(os.path.join(output_dir,'beta_beta_plot_data.tsv.gz')):
  plot_data = pd.DataFrame(
    columns = cnames
  )
  for trait in pipeline_result_files.keys():
    pipeline_res = pd.read_csv(
      os.path.join(pipeline_result_path,pipeline_result_files[trait]),
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
    pipeline_res["Name_pipeline"] = pipeline_res["Name_pipeline"].apply(
      lambda val: val.split(".")[0].strip()
    )
    backman_res = pd.read_csv(
      os.path.join(backman_gwas_path,backman_gwas_files[trait]),
      sep="\t",quoting=csv.QUOTE_NONE
    ).assign(
      Masks = lambda df: df['effect_allele'].map(backman_masks)
    )[["Name","Masks","chromosome",'base_pair_location',"effect_allele","Model","beta","p_value"]].rename(
      {
        "Name":"Name_backman",
        "Model":"Models",
        "base_pair_location":"position_backman",
        "beta":"Beta (Backman et al 2021)",
        "p_value":"Pval (Backman et al 2021)"
      },axis=1
    ).dropna()
    backman_res["Name_backman"] = backman_res["Name_backman"].apply(
      lambda val: re.sub("^.*?\(|\)","",val.split(".")[0])
    )
    shared_data = pd.merge(
      pipeline_res,
      backman_res,
      'inner',
      left_on = ['Masks','chromosome','Models','Name_pipeline'],
      right_on=['Masks','chromosome','Models','Name_backman']
    )
    shared_data['Trait'] = trait
    plot_data = pd.concat(
      [
        plot_data,
        shared_data[cnames]
      ],axis=0
    ).reset_index(drop=True)
    

  with gzip.open(os.path.join(output_dir,'beta_beta_plot_data.tsv.gz'),'wt') as ptr:
    plot_data.to_csv(
      ptr,sep="\t",header=True,index=False,quoting=csv.QUOTE_NONE
    )
else:
  plot_data.read_csv(
      os.path.join(output_dir,'beta_beta_plot_data.tsv.gz'),
      sep="\t",quoting=csv.QUOTE_NONE
    )


row_idx = [x%3 for x in list(range(0,3))]
col_idx = [x%4 for x in list(range(0,4))]
indicies = [(r,c) for r in row_idx for c in col_idx]
fig,ax = plt.subplots(3,4,figsize=(15,15))
for i,t in enumerate(backman_gwas_files.keys()):
  sns.scatterplot(
    plot_data.query(f"Trait == '{t}'"),
    x='Beta (Backman et al 2021)',
    y='Beta (Pipeline results)',ax=ax[indicies[i]]
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
fig.savefig(
  tfile
)
