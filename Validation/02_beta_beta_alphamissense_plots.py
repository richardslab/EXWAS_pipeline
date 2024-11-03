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
output_dir = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/validation_regeneron_10_traits"
# downloaded_data_constants
alphamiss_gwas_path="/scratch/richards/yiheng.chen/project14_ExWAS_AlphaMissense/results/all_regenie_burden_test_res_GWAS_Catalog"
alphamiss_gwas_files = {
  "BMI":"step_2_IRNT_BMI_ExWAS_pLOF_missense_IRNT_BMI.regenie"
}
# own regenie results:
# own regenie results:
pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/pipeline_results/regeneron/Regenie_S2"
pipeline_result_files = {
  "BMI":"8_regenie_S2_OUT_wes_qc_chr*_BMI.regenie.gz"
}


backman_masks = {
  "pLOF_only.singleton" : "M1_LoF.singleton",
  "pLOF_only.01" : "M1_LoF.0.01",
  "pLOF_only.001" : "M1_LoF.0.001",
  'pLOF_and_55missense.singleton' : "M2_LoF_deleterious.singleton",
  'pLOF_and_55missense.01' : "M2_LoF_deleterious.0.01",
  'pLOF_and_55missense.001' : "M2_LoF_deleterious.0.001"
}

cnames = ["Trait","chromosome",'position_backman','position_pipeline',"Masks","Models","Beta (Backman et al 2021)","Beta (Pipeline results)","Pval (Backman et al 2021)","Pval (Pipeline results)","Name_backman","Name_pipeline"]
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
        "Pval":"Pval (Pipeline results)"
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
    os.path.join(alphamiss_gwas_path,alphamiss_gwas_files[trait]),
    sep=" ",quoting=csv.QUOTE_NONE,comment="#"
  ).assign(
    Masks = lambda df: df['ALLELE1'].map(backman_masks),
    p_value = lambda df: df['LOG10P'].apply(lambda val: 10**(-1 * val)),
    Models = lambda df: df['TEST'].apply(lambda val: f"{val}-WGR-LR")
  )[["ID","Masks","CHROM",'GENPOS',"ALLELE1","Models","BETA","p_value"]].rename(
    {
      "ID":"Name_backman",
      "CHROM":"chromosome",
      "GENPOS":"position_backman",
      "BETA":"Beta (Backman et al 2021)",
      "p_value":"Pval (Backman et al 2021)"
    },axis=1
  ).dropna()
  backman_res["Name_backman"] = backman_res["Name_backman"].apply(
    lambda val: re.sub("\..*$","",val)
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



row_idx = [x%3 for x in list(range(0,3))]
col_idx = [x%4 for x in list(range(0,4))]
indicies = [(r,c) for r in row_idx for c in col_idx]
fig,ax = plt.subplots(3,4,figsize=(15,15))
for i,t in enumerate(alphamiss_gwas_files.keys()):
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
fig.savefig(tfile)
