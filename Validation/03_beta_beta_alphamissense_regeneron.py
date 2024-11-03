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
backman_gwas_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/data/backman_exwas_results"
backman_gwas_files = {
  "BMI":"GCST90082670_buildGRCh38.tsv.gz"
}

backman_masks = {
  "M1.singleton" : "pLOF_only.singleton",
  "M1.01" : "pLOF_only.01",
  "M1.001" : "pLOF_only.001",
  'M3.singleton' : "pLOF_and_55missense.singleton",
  'M3.01' : "pLOF_and_55missense.01",
  'M3.001' : "pLOF_and_55missense.001"
}

cnames = ["Trait","chromosome",'position_backman','position_alphamiss',"Masks","Models","Beta (Backman et al 2021)","Beta (alphamiss)","Pval (Backman et al 2021)","Pval (alphamiss)","Name_backman","Name_alphamiss"]
# Compare results

plot_data = pd.DataFrame(
  columns = cnames
)
for trait in pipeline_result_files.keys():
  alphamiss_res = pd.read_csv(
    os.path.join(alphamiss_gwas_path,alphamiss_gwas_files[trait]),
    sep=" ",quoting=csv.QUOTE_NONE,comment="#"
  ).assign(
    p_value = lambda df: df['LOG10P'].apply(lambda val: 10**(-1 * val)),
    Models = lambda df: df['TEST'].apply(lambda val: f"{val}-WGR-LR")
  )[["ID","CHROM",'GENPOS',"ALLELE1","Models","BETA","p_value"]].rename(
    {
      "ID":"Name_alphamiss",
      "CHROM":"chromosome",
      "ALLELE1":"Masks",
      "GENPOS":"position_alphamiss",
      "BETA":"Beta (alphamiss)",
      "p_value":"Pval (alphamiss)"
    },axis=1
  ).dropna()
  alphamiss_res["Name_alphamiss"] = alphamiss_res["Name_alphamiss"].apply(
    lambda val: re.sub("\..*$","",val)
  )
  # backman
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
  # merge
  shared_data = pd.merge(
    alphamiss_res,
    backman_res,
    'inner',
    left_on = ['Masks','chromosome','Models','Name_alphamiss'],
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
    plot_data.query(f"Trait == '{t}' & chromosome == 2"),
    x='Beta (Backman et al 2021)',
    y='Beta (alphamiss)',ax=ax[indicies[i]]
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
fig.supylabel("alphamiss")
fig.savefig(tfile)


chr21 = plot_data.query("chromosome == 21")
dup_genes = chr21.value_counts(['Name_backman']).to_frame().set_axis(['count'],axis=1)
chr21.query("Name_backman == 'ENSG00000141956'")