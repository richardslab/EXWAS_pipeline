"""
  Compare the effect sizes my runs and Alphamissense
"""
#%%

import os,shutil,yaml,pickle,re,glob,csv,gzip,math
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
ffile_bb_plof = os.path.join(outdir,'beta_beta_alphamissense_plof_binary.png')
ffile_pp_plof = os.path.join(outdir,'pval_pval_alphamissense_plof_binary.png')

ffile_bb_plof_5in5 = os.path.join(outdir,'beta_beta_alphamissense_plof_or_d5in5_binary.png')
ffile_pp_plof_5in5 = os.path.join(outdir,'pval_pval_alphamissense_plof_or_5in5_binary.png')

# downloaded_data_constants
alphamiss_gwas_path="/scratch/richards/yiheng.chen/project14_ExWAS_AlphaMissense/results/all_regenie_burden_test_res_GWAS_Catalog"
alphamiss_gwas_files = {
  "hypertension":"step_2_hypertension_ExWAS_pLOF_missense_hypertension.regenie",
  "Hypercholesterolemia":"step_2_Hypercholesterolemia_ExWAS_pLOF_missense_Hypercholesterolemia.regenie",
  "Cataract":"step_2_Cataract_ExWAS_pLOF_missense_Cataract.regenie",
  "T2D":"step_2_T2D_ExWAS_pLOF_missense_T2D.regenie",
  "Hypothyroidism":"step_2_Hypothyroidism_ExWAS_pLOF_missense_Hypothyroidism.regenie",
  "Acute_renal_failure":"step_2_Acute_renal_failure_ExWAS_pLOF_missense_Acute_renal_failure.regenie",
  "Atrial_fibrillation":"step_2_Atrial_fibrillation_ExWAS_pLOF_missense_Atrial_fibrillation.regenie",
  "Osteoarthritis_localized":"step_2_Osteoarthritis_localized_ExWAS_pLOF_missense_Osteoarthritis_localized.regenie",
  "Diaphragmatic_hernia":"step_2_Diaphragmatic_hernia_ExWAS_pLOF_missense_Diaphragmatic_hernia.regenie",
  "MDD":"step_2_Major_depressive_disorder_ExWAS_pLOF_missense_Major_depressive_disorder.regenie"
}
# own regenie results:
pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/alphamiss_exact_binary/alphamiss_plof_5in5/Regenie_S2"
pipeline_result_files = {
  "hypertension":"8_regenie_S2_OUT_ukb_merged_1-22_hypertension.regenie.gz",
  "Hypercholesterolemia":"8_regenie_S2_OUT_ukb_merged_1-22_Hypercholesterolemia.regenie.gz",
  "Cataract":"8_regenie_S2_OUT_ukb_merged_1-22_Cataract.regenie.gz",
  "T2D":"8_regenie_S2_OUT_ukb_merged_1-22_T2D.regenie.gz",
  "Hypothyroidism":"8_regenie_S2_OUT_ukb_merged_1-22_Hypothyroidism.regenie.gz",
  "Acute_renal_failure":"8_regenie_S2_OUT_ukb_merged_1-22_Acute_renal_failure.regenie.gz",
  "Atrial_fibrillation":"8_regenie_S2_OUT_ukb_merged_1-22_Atrial_fibrillation.regenie.gz",
  "Osteoarthritis_localized":"8_regenie_S2_OUT_ukb_merged_1-22_Osteoarthritis_localized.regenie.gz",
  "Diaphragmatic_hernia":"8_regenie_S2_OUT_ukb_merged_1-22_Diaphragmatic_hernia.regenie.gz",
  "MDD":"8_regenie_S2_OUT_ukb_merged_1-22_Major_depressive_disorder.regenie.gz"
}

Alphamissense_masks = {
  "pLOF_only.singleton": "M1_LoF.singleton",
  "pLOF_only.0.01": "M1_LoF.0.01",
  "pLOF_only.0.001": "M1_LoF.0.001",
  "pLOF_and_55missense.singleton": 'M2_LoF_or_deleterious.singleton',
  "pLOF_and_55missense.0.01": 'M2_LoF_or_deleterious.0.01',
  "pLOF_and_55missense.0.001": 'M2_LoF_or_deleterious.0.001'
}




cnames = [
  "Trait","chromosome",'position_Alphamissense','position_pipeline',"Masks","Models",
  "Beta (Chen et al 2024)","Beta (Pipeline results)",
  "Pval (Chen et al 2024)","Pval (Pipeline results)",
  "LOG10P (Chen et al 2024)","LOG10P (Pipeline results)",
  "SE (Chen et al 2024)","SE (Pipeline results)",
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

def _parse_info_beta(val):
  res = re.findall("REGENIE_BETA=.*?;",val)
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
      ),
      beta = lambda df: df['Info'].apply(
        lambda val: _parse_info_beta(val)
      )
    )[['Name',"Alt","Model","Chr",'Pos',"Ref","Pval","beta",'LOG10P','SE']].rename(
      {
        "Name":"Name_pipeline",
        "Alt":"Masks",
        "Chr":"chromosome",
        "Pos":"position_pipeline",
        "Model":"Models",
        "beta":"Beta (Pipeline results)",
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
    Models = lambda df: df['TEST'].apply(lambda val: f"{val}-WGR-FIRTH")
  )[["ID","Masks","CHROM",'GENPOS',"ALLELE1","Models","BETA","p_value","SE",'LOG10P']].rename(
    {
      "ID":"Name_Alphamissense",
      "CHROM":"chromosome",
      "GENPOS":"position_Alphamissense",
      "BETA":"Beta (Chen et al 2024)",
      "p_value":"Pval (Chen et al 2024)",
      "SE":"SE (Chen et al 2024)",
      "LOG10P":"LOG10P (Chen et al 2024)"
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

def make_fig(plot_data,x,y,title,beta=False):
  nrow=3
  ncol=4
  row_idx = [x%nrow for x in list(range(0,nrow))]
  col_idx = [x%ncol for x in list(range(0,ncol))]
  indicies = [(r,c) for r in row_idx for c in col_idx]
  fig,ax = plt.subplots(nrow,ncol,figsize=(33,20))
  for i,t in enumerate(alphamiss_gwas_files.keys()):
    trait_plot = plot_data.query(f"Trait == '{t}'").drop(['Masks','Models'],axis=1)
    trait_plot = trait_plot.drop_duplicates()
    if beta:
      min_scale = math.floor(min(trait_plot[x].min(),trait_plot[y].min()))-5
    else:
      min_scale = max(0,math.floor(min(trait_plot[x].min(),trait_plot[y].min()))-5)
    max_scale = math.ceil(max(trait_plot[x].max(),trait_plot[y].max())) + 5
    scale = [math.ceil(x) - math.ceil(x)%5 for x in np.linspace(min_scale,max_scale,5)]
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
        "title":f"{t} "+rf"$R^2$: {np.round(r2.correlation,2)}"
      }
    )
    ax[indicies[i]].set_xticks(scale)
    ax[indicies[i]].set_yticks(scale)
    ax[indicies[i]].set_xlim(min_scale,max_scale)
    ax[indicies[i]].set_ylim(min_scale,max_scale)
  fig.supxlabel("Chen et al 2024")
  fig.supylabel("Pipeline results")
  fig.suptitle(title)
  return(fig)



# Pval plof
fig_plof = make_fig(
  plot_data=plot_data.query("Masks.str.startswith('M1_LoF')"),
  x = 'LOG10P (Chen et al 2024)',
  y = "LOG10P (Pipeline results)",
  title = r"$log_{10}(P-value)$ vs $log_{10}(P-value)$"+"\nChen et al 2024 (pLoF)"
)
fig_plof.savefig(tfile)
fig_plof.savefig(f"{ffile_pp_plof}")
plt.close(fig_plof)

# Pval plof or 5in5
fig_plof = make_fig(
  plot_data=plot_data.query("Masks.str.startswith('M2_')"),
  x = 'LOG10P (Chen et al 2024)',
  y = "LOG10P (Pipeline results)",
  title = r"$log_{10}(P-value)$ vs $log_{10}(P-value)$"+"\nChen et al 2024 (pLoF or deleterious 5 in 5)"
)
fig_plof.savefig(tfile)
fig_plof.savefig(f"{ffile_pp_plof_5in5}")
plt.close(fig_plof)


