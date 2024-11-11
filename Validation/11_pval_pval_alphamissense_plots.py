"""
  Compare the effect sizes my runs and Alphamissense
"""
#%%

import os,shutil,yaml,pickle,re,glob,csv,gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import tempfile
from scipy import stats

tempdir = tempfile.TemporaryDirectory()
tfile = os.path.join(tempdir.name,'temp.png')
outdir = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/validation_regeneron_10_traits"
ffile = os.path.join(outdir,'pval_pval_alphamissense_plof_or_5in5.pdf')
plt.rcParams.update(
  {
    "font.size":25
  }
)


#
output_dir = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/validation_regeneron_10_traits" 
# downloaded_data_constants
alphamiss_gwas_path="/scratch/richards/yiheng.chen/project14_ExWAS_AlphaMissense/results/all_regenie_burden_test_res_GWAS_Catalog"
alphamiss_gwas_files = {
  "ZBMD":"step_2_ZBMD_ExWAS_pLOF_missense_ZBMD.regenie",
  "dBilirubin":"step_2_IRNT_biliru_ExWAS_pLOF_missense_IRNT_biliru.regenie",
  "Calcium":"step_2_IRNT_Ca_ExWAS_pLOF_missense_IRNT_Ca.regenie",
  "DBP":"step_2_IRNT_DBP_ExWAS_pLOF_missense_IRNT_DBP.regenie",
  "Glucose":"step_2_IRNT_glu_ExWAS_pLOF_missense_IRNT_glu.regenie",
  "Height":"step_2_IRNT_height_ExWAS_pLOF_missense_IRNT_height.regenie",
  "LDL":"step_2_IRNT_LDL_ExWAS_pLOF_missense_IRNT_LDL.regenie",
  "RBC":"step_2_IRNT_RBC_ExWAS_pLOF_missense_IRNT_RBC.regenie",
  "SBP":"step_2_IRNT_SBP_ExWAS_pLOF_missense_IRNT_SBP.regenie",
  "TG":"step_2_IRNT_TG_ExWAS_pLOF_missense_IRNT_TG.regenie",
  "WHR":"step_2_IRNT_WHR_ExWAS_pLOF_missense_IRNT_WHR.regenie"
}
# own regenie results:
# own regenie results:
pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/pipeline_results_plof_or_5in5/regeneron/Regenie_S2"
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


Alphamissense_masks = {
  "pLOF_only.singleton" : "M1_LoF.singleton",
  "pLOF_only.01" : "M1_LoF.0.01",
  "pLOF_only.001" : "M1_LoF.0.001",
  'pLOF_and_55missense.singleton' : "M2_LoF_deleterious.singleton",
  'pLOF_and_55missense.01' : "M2_LoF_deleterious.0.01",
  'pLOF_and_55missense.001' : "M2_LoF_deleterious.0.001"
}

cnames = ["Trait","chromosome",'position_Alphamissense','position_pipeline',"Masks","Models","Beta (Alphamissense)","Beta (Pipeline results)","Pval (Alphamissense)","Pval (Pipeline results)","Name_Alphamissense","Name_pipeline","Log10P (Pipeline results)","Log10P (Alphamissense)"]
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
      sep=" ",quoting=csv.QUOTE_NONE,comment='#'
    ).assign(
      Models = lambda df: df['TEST'].apply(lambda val: f"{val}-WGR-LR"),
      Pval = lambda df: df['LOG10P'].apply(lambda val: 10**(-1 * val))
    )[['ID',"ALLELE1","Models","CHROM",'GENPOS',"ALLELE0","LOG10P","Pval","BETA"]].rename(
      {
        "ID":"Name_pipeline",
        "ALLELE1":"Masks",
        "CHROM":"chromosome",
        "GENPOS":"position_pipeline",
        "BETA":"Beta (Pipeline results)",
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
  Alphamissense_res = pd.read_csv(
    os.path.join(alphamiss_gwas_path,alphamiss_gwas_files[trait]),
    sep=" ",quoting=csv.QUOTE_NONE,comment="#"
  ).assign(
    Masks = lambda df: df['ALLELE1'].map(Alphamissense_masks),
    p_value = lambda df: df['LOG10P'].apply(lambda val: 10**(-1 * val)),
    Models = lambda df: df['TEST'].apply(lambda val: f"{val}-WGR-LR")
  )[["ID","Masks","CHROM",'GENPOS',"ALLELE1","Models","BETA","p_value",'LOG10P']].rename(
    {
      "ID":"Name_Alphamissense",
      "CHROM":"chromosome",
      "GENPOS":"position_Alphamissense",
      "BETA":"Beta (Alphamissense)",
      "p_value":"Pval (Alphamissense)",
      "LOG10P":"Log10P (Alphamissense)",
    },axis=1
  ).dropna()
  Alphamissense_res["Name_Alphamissense"] = Alphamissense_res["Name_Alphamissense"].apply(
    lambda val: re.sub("\..*$","",val)
  )
  shared_data = pd.merge(
    pipeline_results,
    Alphamissense_res,
    'inner',
    left_on = ['Masks','chromosome','Models','Name_pipeline'],
    right_on=['Masks','chromosome','Models','Name_Alphamissense']
  )
  shared_data['Trait'] = trait
  plot_data = pd.concat(
    [
      plot_data,
      shared_data[cnames]
    ],axis=0
  ).reset_index(drop=True)



row_idx = [x%3 for x in list(range(0,3))]
col_idx = [x%4 for x in list(range(0,4))]
indicies = [(r,c) for r in row_idx for c in col_idx]
fig,ax = plt.subplots(3,4,figsize=(25,20))
for i,t in enumerate(alphamiss_gwas_files.keys()):
  trait_plot = plot_data.query(f"Trait == '{t}'")
  r2 = stats.pearsonr(
    trait_plot['Log10P (Alphamissense)'],
    trait_plot['Log10P (Pipeline results)']
  )
  sns.scatterplot(
    trait_plot,
    x='Log10P (Alphamissense)',
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

fig.supxlabel(r"$log_{10}(p-value)$ Alphamissense")
fig.supylabel(r"$log_{10}(p-value)$ Pipeline results")
fig.savefig(tfile)


fig.savefig(ffile)