"""
  Compare the effect sizes between runs.

  The 2 ways to specify masks are the same
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

ffile_pp_plof = os.path.join(outdir,'pval_pval_self.png')

# downloaded_data_constants
pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/alphamiss_exact/alphamiss_plof/Regenie_S2"
pipeline_result_files = {
  "ZBMD":"8_regenie_S2_OUT_ukb_merged_1-22_ZBMD.regenie.gz",
  "dBilirubin":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_biliru.regenie.gz",
  "Calcium":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_Ca.regenie.gz",
  "DBP":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_DBP.regenie.gz",
  "Glucose":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_glu.regenie.gz",
  "Height":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_height.regenie.gz",
  "LDL":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_LDL.regenie.gz",
  "RBC":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_RBC.regenie.gz",
  "SBP":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_SBP.regenie.gz",
  "TG":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_TG.regenie.gz",
  "WHR":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_WHR.regenie.gz"
}


# own regenie results:
pipeline_result_path_2="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/alphamiss_exact/alphamiss_plof_5in5/Regenie_S2"
pipeline_result_files_2 = {
  "ZBMD":"8_regenie_S2_OUT_ukb_merged_1-22_ZBMD.regenie.gz",
  "dBilirubin":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_biliru.regenie.gz",
  "Calcium":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_Ca.regenie.gz",
  "DBP":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_DBP.regenie.gz",
  "Glucose":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_glu.regenie.gz",
  "Height":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_height.regenie.gz",
  "LDL":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_LDL.regenie.gz",
  "RBC":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_RBC.regenie.gz",
  "SBP":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_SBP.regenie.gz",
  "TG":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_TG.regenie.gz",
  "WHR":"8_regenie_S2_OUT_ukb_merged_1-22_IRNT_WHR.regenie.gz"
}


cnames = [
  "Trait","chromosome",'position_pipeline','position_pipeline_1',"Masks","Models",
  "Beta (Pipeline results)","Beta (Pipeline results 1)",
  "Pval (Pipeline results)","Pval (Pipeline results 1)",
  "LOG10P (Pipeline results)","LOG10P (Pipeline results 1)",
  "SE (Pipeline results)","SE (Pipeline results 1)",
  "Name_pipeline","Name_pipeline_1"
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
  all_pipeline_results_2 = glob.glob(os.path.join(pipeline_result_path_2,pipeline_result_files_2[trait]))
  pipeline_results_2 = pd.DataFrame()
  for each_f in all_pipeline_results_2:
    each_res = pd.read_csv(
      os.path.join(pipeline_result_path_2,each_f),
      sep="\t",quoting=csv.QUOTE_NONE,comment='#'
    ).query("~Effect.isna()").assign(
      LOG10P = lambda df: df['Pval'].apply(lambda val: -1 * np.log10(val)),
      SE = lambda df: df['Info'].apply(
        lambda val: _parse_info(val)
      )
    )[['Name',"Alt","Model","Chr",'Pos',"Ref","Pval","Effect",'LOG10P','SE']].rename(
      {
        "Name":"Name_pipeline_1",
        "Alt":"Masks",
        "Chr":"chromosome",
        "Pos":"position_pipeline_1",
        "Model":"Models",
        "Effect":"Beta (Pipeline results 1)",
        "Pval":"Pval (Pipeline results 1)",
        "LOG10P":"LOG10P (Pipeline results 1)",
        "SE":"SE (Pipeline results 1)"
      },axis=1
    )
    each_res["Name_pipeline_1"] = each_res["Name_pipeline_1"].apply(
      lambda val: val.split(".")[0].strip()
    )
    pipeline_results_2 = pd.concat(
      [
        pipeline_results_2,
        each_res
      ],axis=0
    ).reset_index(drop=True)
  shared_data = pd.merge(
    pipeline_results,
    pipeline_results_2,
    'inner',
    left_on = ['Masks','chromosome','Models','Name_pipeline'],
    right_on=['Masks','chromosome','Models','Name_pipeline_1']
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
  fig,ax = plt.subplots(nrow,ncol,figsize=(30,30))
  for i,t in enumerate(pipeline_result_files.keys()):
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
        "title":f"{t}\n"+rf"$R^2$: {np.round(r2.correlation,3)}"
      }
    )
    ax[indicies[i]].set_xticks(scale)
    ax[indicies[i]].set_yticks(scale)
    ax[indicies[i]].set_xlim(min_scale,max_scale)
    ax[indicies[i]].set_ylim(min_scale,max_scale)
  fig.supxlabel("Pipeline results 1")
  fig.supylabel("Pipeline results")
  fig.suptitle(title)
  return(fig)

fig_plof = make_fig(
  plot_data=plot_data.query("Masks.str.startswith('M1_LoF')"),
  x = 'Beta (Pipeline results 1)',
  y = "Beta (Pipeline results)",
  title = r"$\beta$ vs $\beta$"+"\Pipeline results 1 (pLoF)",
  beta=True
)
fig_plof.savefig(f"{tfile}")
fig_plof.savefig(f"{ffile_pp_plof}")
