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
import statsmodels.api as sm
from scipy import stats

tempdir = tempfile.TemporaryDirectory()
tfile = os.path.join(tempdir.name,'temp.png')
plt.rcParams.update(
  {
    "font.size":25
  }
)
outdir = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/validation_regeneron_10_traits"
ffile = os.path.join(outdir,'pval_pval_plots_regeneron_plof_5in5.pdf')


# downloaded_data_constants
backman_gwas_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/data/backman_exwas_results"
backman_gwas_files = {
  "ZBMD":"GCST90083038_buildGRCh38.tsv.gz",
  "dBilirubin":"GCST90083007_buildGRCh38.tsv.gz",
  "Calcium":"GCST90083009_buildGRCh38.tsv.gz",
  "DBP":"GCST90083131_buildGRCh38.tsv.gz",
  "Glucose":"GCST90083015_buildGRCh38.tsv.gz",
  "Height":"GCST90083263_buildGRCh38.tsv.gz",
  "LDL":"GCST90083019_buildGRCh38.tsv.gz",
  "RBC":"GCST90082965_buildGRCh38.tsv.gz",
  "SBP":"GCST90083132_buildGRCh38.tsv.gz",
  "TG":"GCST90083030_buildGRCh38.tsv.gz"
}
# own regenie results:
pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/PAST_pipeline_results/validation_10_traits_plof_or_5in5/regeneron/Regenie_S2"
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
  "TG":"8_regenie_S2_OUT_wes_qc_chr*_IRNT_TG.regenie.gz"
}

backman_masks = {
  "M1.singleton" : "M1_LoF.singleton",
  "M1.01" : "M1_LoF.0.01",
  "M1.001" : "M1_LoF.0.001",
  "M1.0001" : "M1_LoF.0.0001",
  'M3.singleton' : "M2_LoF_deleterious.singleton",
  'M3.01' : "M2_LoF_deleterious.0.01",
  'M3.001' : "M2_LoF_deleterious.0.001",
  'M3.0001' : "M2_LoF_deleterious.0.0001"
}

cnames = ["Trait","chromosome",'position_backman','position_pipeline',"Masks","Models","Beta (Backman et al 2021)","Beta (Pipeline results)","Pval (Backman et al 2021)","Pval (Pipeline results)","Name_backman","Name_pipeline","Log10P (Backman et al 2021)","Log10P (Pipeline results)"]
# Compare results
plot_data = pd.DataFrame(
    columns = cnames
  )
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
    pipeline_results = pd.concat(
      [
        pipeline_results,
        each_res
      ],axis=0
    ).reset_index(drop=True)
  backman_res = pd.read_csv(
    os.path.join(backman_gwas_path,backman_gwas_files[trait]),
    sep="\t",quoting=csv.QUOTE_NONE
  ).assign(
    Masks = lambda df: df['effect_allele'].map(backman_masks)
  )[["Name","Masks","chromosome",'base_pair_location',"effect_allele","Model","beta","p_value"]].assign(
    log10p = lambda df: -1 * np.log10(df['p_value'])
  ).rename(
    {
      "Name":"Name_backman",
      "Model":"Models",
      "base_pair_location":"position_backman",
      "beta":"Beta (Backman et al 2021)",
      "p_value":"Pval (Backman et al 2021)",
      "log10p": "Log10P (Backman et al 2021)"
    },axis=1
  ).dropna()
  backman_res["Name_backman"] = backman_res["Name_backman"].apply(
    lambda val: re.sub("^.*?\(|\)","",val.split(".")[0])
  )
  shared_data = pd.merge(
    pipeline_results,
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


row_idx = [x%2 for x in list(range(0,2))]
col_idx = [x%5 for x in list(range(0,5))]
indicies = [(r,c) for r in row_idx for c in col_idx]
fig,ax = plt.subplots(2,5,figsize=(30,15))
for i,t in enumerate(backman_gwas_files.keys()):
  trait_plot = plot_data.query(f"Trait == '{t}' & `Pval (Backman et al 2021)` < 0.05")
  # get correlation
  r2 = stats.pearsonr(
    trait_plot['Log10P (Backman et al 2021)'],
    trait_plot['Log10P (Pipeline results)']
  )
  sns.scatterplot(
    trait_plot,
    x='Log10P (Backman et al 2021)',
    y='Log10P (Pipeline results)',ax=ax[indicies[i]]
  )
  ax[indicies[i]].axline(
    (0,0),slope = 1,linestyle='--'
  )
  ax[indicies[i]].update(
    {
      "xlabel":"",
      "ylabel":"",
      "title":f"{t} R2: {np.round(r2.correlation,2)}"
    }
  )

fig.supxlabel(r"$-log_{10}(p-value)$ Backman et al 2021")
fig.supylabel(r"$-log_{10}(p-value)$ Pipeline results")
fig.savefig(tfile)

fig.savefig(ffile)



