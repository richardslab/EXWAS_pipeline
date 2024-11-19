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

plt.rcParams.update(
  {
    "font.size":25
  }
)
outdir = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/figures"
ffile_bb = os.path.join(outdir,'beta_beta_alphamissense_plof.png')
ffile_pp = os.path.join(outdir,'pval_pval_alphamissense_plof.png')

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
pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/plof/regeneron/Regenie_S2"
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
  "pLOF_only.0.01" : "M1_LoF.0.01",
  "pLOF_only.0.001" : "M1_LoF.0.001",
  'pLOF_and_55missense.singleton' : "M4_LoF_or_deleterious_5_of_5.singleton",
  'pLOF_and_55missense.0.01' : "M4_LoF_or_deleterious_5_of_5.0.01",
  'pLOF_and_55missense.0.001' : "M4_LoF_or_deleterious_5_of_5.0.001"
}



cnames = [
  "Trait","chromosome",'position_Alphamissense','position_pipeline',"Masks","Models",
  "Beta (Alphamissense)","Beta (Pipeline results)",
  "Pval (Alphamissense)","Pval (Pipeline results)",
  "LOG10P (Alphamissense)","LOG10P (Pipeline results)",
  "SE (Alphamissense)","SE (Pipeline results)",
  "Name_Alphamissense","Name_pipeline"
]

# Compare results
def _parse_info(val):
  res = re.findall("REGENIE_SE=.*?;",val)
  if len(res) == 1:
    res = float(res[0].split("=")[1].strip(";"))
  else:
    res = None
  return res

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
  )[["ID","Masks","CHROM",'GENPOS',"ALLELE1","Models","BETA","p_value","SE",'LOG10P']].rename(
    {
      "ID":"Name_Alphamissense",
      "CHROM":"chromosome",
      "GENPOS":"position_Alphamissense",
      "BETA":"Beta (Alphamissense)",
      "p_value":"Pval (Alphamissense)",
      "SE":"SE (Alphamissense)",
      "LOG10P":"LOG10P (Alphamissense)"
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

def make_fig(filtered,x,y,title):
  row_idx = [x%3 for x in list(range(0,3))]
  col_idx = [x%4 for x in list(range(0,4))]
  indicies = [(r,c) for r in row_idx for c in col_idx]
  fig,ax = plt.subplots(3,4,figsize=(30,20))
  for i,t in enumerate(alphamiss_gwas_files.keys()):
    if filtered:
      trait_plot = plot_data.query(f"Trait == '{t}' & `Pval (Alphamissense)` < 0.05 & `Pval (Pipeline results)` < 0.05")
    else:
      trait_plot = plot_data.query(f"Trait == '{t}'")
    # get correlation
    r2 = stats.pearsonr(
      trait_plot[x],
      trait_plot[y]
    )
    sns.scatterplot(
      trait_plot,
      x=x,
      y=y,ax=ax[indicies[i]]
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
  fig.supxlabel("Alphamissense")
  fig.supylabel("Pipeline results")
  fig.suptitle(title)
  return(fig)

fig = make_fig(
  filtered = False,
  x = 'Beta (Alphamissense)',
  y = "Beta (Pipeline results)",
  title = r"$\beta$ - $\beta$ plots Alphamissense (pLoF)"
)
# fig.savefig(tfile)
fig.savefig(f"{ffile_bb}")
plt.close(fig)

fig = make_fig(
  filtered = True,
  x = 'Beta (Alphamissense)',
  y = "Beta (Pipeline results)",
  title = r"$\beta$ - $\beta$ plots Alphamissense (pLoF)"
)
# fig.savefig(tfile)
fig.savefig(f"{ffile_bb}_filtered.png")
plt.close(fig)

# Pval
fig = make_fig(
  filtered = False,
  x = 'LOG10P (Alphamissense)',
  y = "LOG10P (Pipeline results)",
  title = r"$log_{10}(P-value)$ - $log_{10}(P-value)$ plots Alphamissense (pLoF)"
)
fig.savefig(tfile)
fig.savefig(f"{ffile_pp}")
plt.close(fig)

fig = make_fig(
  filtered = True,
  x = 'LOG10P (Alphamissense)',
  y = "LOG10P (Pipeline results)",
  title = r"$log_{10}(P-value)$ - $log_{10}(P-value)$ plots Alphamissense (pLoF)"
)
fig.savefig(tfile)
fig.savefig(f"{ffile_pp}_filtered.png")
plt.close(fig)


