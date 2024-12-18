"""
Obtain the list of hits found by Regeneron and by the pipelien based on their critiera
"""

import os,shutil,yaml,re,csv,scipy,glob
import pandas as pd
import numpy as np

outdir = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/tables"

#%% pipeline backman results
pipeline_backman_path = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/regeneron_exact/regeneron/Regenie_S2"
pipeline_backman_files = {
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
def _parse_info(val):
  res = re.findall("REGENIE_SE=.*?;",val)
  if len(res) == 1:
    res = float(res[0].split("=")[1].strip(";"))
  else:
    res = None
  return res

all_pipeline_backman_res = pd.DataFrame()
for pheno,file in pipeline_backman_files.items():
  all_files = glob.glob(os.path.join(pipeline_backman_path,file))
  for each_f in all_files:
    each_res = pd.read_csv(
      each_f,
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
    each_res['Trait'] = pheno
    all_pipeline_backman_res = pd.concat(
      [
        all_pipeline_backman_res,
        each_res
      ],axis=0
    ).reset_index(drop=True)

# pipeline chen results
pipeline_alpha_continuous_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/alphamiss_exact/alphamiss_plof_5in5/Regenie_S2"
pipeline_alpha_continuous_files = {
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
all_pipeline_alpha_res = pd.DataFrame()
for pheno,file in pipeline_alpha_continuous_files.items():
  all_files = glob.glob(os.path.join(pipeline_alpha_continuous_path,file))
  for each_f in all_files:
    each_res = pd.read_csv(
      os.path.join(pipeline_alpha_continuous_path,each_f),
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
    each_res['Trait'] = pheno
    all_pipeline_alpha_res = pd.concat(
      [
        all_pipeline_alpha_res,
        each_res
      ],axis=0
    ).reset_index(drop=True)

pipeline_alpha_binary_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/alphamiss_exact_binary/alphamiss_plof_5in5/Regenie_S2"
pipeline_alpha_binary_files = {
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

def _parse_info_beta(val):
  res = re.findall("REGENIE_BETA=.*?;",val)
  if len(res) == 1:
    res = float(res[0].split("=")[1].strip(";"))
  else:
    res = None
  return res

for pheno,file in pipeline_alpha_binary_files.items():
  all_files = glob.glob(os.path.join(pipeline_alpha_binary_path,file))
  for each_f in all_files:
    each_res = pd.read_csv(
      os.path.join(pipeline_alpha_binary_path,each_f),
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
    each_res['Trait'] = pheno
    all_pipeline_alpha_res = pd.concat(
      [
        all_pipeline_alpha_res,
        each_res
      ],axis=0
    ).reset_index(drop=True)

# save

all_pipeline_backman_res.to_csv(
  os.path.join(
    outdir,'all_pipeline_backman_res.tsv.gz'
  ),sep="\t",quoting=csv.QUOTE_NONE,index=False
)

all_pipeline_alpha_res.to_csv(
  os.path.join(
    outdir,'all_pipeline_alpha_res.tsv.gz'
  ),sep="\t",quoting=csv.QUOTE_NONE,index=False
)