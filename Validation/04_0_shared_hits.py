"""
Obtain the list of hits found by Regeneron and by the pipelien based on their critiera
"""

import os,shutil,yaml,re,csv,scipy,glob
import pandas as pd
import numpy as np

outdir = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/tables"

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

pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/plof_or_5in5/regeneron/Regenie_S2"
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

# obtain all pipeline results
def _parse_info(val):
  res = re.findall("REGENIE_SE=.*?;",val)
  if len(res) == 1:
    res = float(res[0].split("=")[1].strip(";"))
  else:
    res = None
  return res

pipeline_results = pd.DataFrame()

for trait in pipeline_result_files.keys():
  all_pipeline_results = glob.glob(os.path.join(pipeline_result_path,pipeline_result_files[trait]))
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
    each_res['Trait'] = trait
    pipeline_results = pd.concat(
      [
        pipeline_results,
        each_res
      ],axis=0
    ).reset_index(drop=True)

# backman results
backman_masks = {
  "M1.singleton" : "M1_LoF.singleton",
  "M1.0001" : "M1_LoF.1e-05",
  "M1.001" : "M1_LoF.0.0001",
  "M1.01" : "M1_LoF",
  "M1.1" : "M1_LoF.0.01",
  'M3.singleton' : "M4_LoF_or_deleterious_5_of_5.singleton",
  'M3.0001' : "M4_LoF_or_deleterious_5_of_5.1e-05",
  'M3.001' : "M4_LoF_or_deleterious_5_of_5.0.0001",
  'M3.01' : "M4_LoF_or_deleterious_5_of_5.0.001",
  'M3.1' : "M4_LoF_or_deleterious_5_of_5.0.01"
}

backman_results = pd.DataFrame()

for trait in pipeline_result_files.keys():
  backman_res = pd.read_csv(
    os.path.join(backman_gwas_path,backman_gwas_files[trait]),
    sep="\t",quoting=csv.QUOTE_NONE
  ).assign(
    Masks = lambda df: df['effect_allele'].map(backman_masks),
    LOG10P = lambda df: df['p_value'].apply(lambda val: -1 * np.log10(val))
  )[["Name","Masks","chromosome",'base_pair_location',"effect_allele","Model","beta","p_value",'LOG10P','standard_error']].rename(
    {
      "Name":"Name_backman",
      "Model":"Models",
      "base_pair_location":"position_backman",
      "beta":"Beta (Backman et al 2021)",
      "p_value":"Pval (Backman et al 2021)",
      "LOG10P":"LOG10P (Backman et al 2021)",
      "standard_error":"SE (Backman et al 2021)"
    },axis=1
  ).dropna()
  backman_res["Name_backman"] = backman_res["Name_backman"].apply(
    lambda val: re.sub("^.*?\(|\)","",val.split(".")[0])
  )
  backman_res = backman_res.query("Masks.str.contains('M4_LoF_or_deleterious_5_of_5')")
  backman_res['Trait'] = trait
  backman_results = pd.concat(
    [backman_results,backman_res],axis=0
  ).reset_index(drop=True)

# alphamissense results
Alphamissense_masks = {
  "pLOF_only.singleton" : "M1_LoF.singleton",
  "pLOF_only.0.01" : "M1_LoF.0.01",
  "pLOF_only.0.001" : "M1_LoF.0.001",
  'pLOF_and_55missense.singleton' : "M4_LoF_or_deleterious_5_of_5.singleton",
  'pLOF_and_55missense.0.01' : "M4_LoF_or_deleterious_5_of_5.0.01",
  'pLOF_and_55missense.0.001' : "M4_LoF_or_deleterious_5_of_5.0.001"
}

alphamissense_results = pd.DataFrame()

for trait in pipeline_result_files.keys():
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
  Alphamissense_res["Trait"] = trait
  Alphamissense_res = Alphamissense_res.query("Masks.str.contains('M4_LoF_or_deleterious_5_of_5')")
  alphamissense_results = pd.concat(
    [
      alphamissense_results,
      Alphamissense_res 
    ],axis=0
  ).reset_index(drop=True)

pipeline_alpha = pd.merge(
  pipeline_results,alphamissense_results,
  'inner',
  left_on=['Name_pipeline','Masks',"Models","chromosome",'Trait'],
  right_on=['Name_Alphamissense','Masks',"Models","chromosome",'Trait']
)
pipeline_regeneron = pd.merge(
  pipeline_results,backman_results,
  'inner',
  left_on=['Name_pipeline','Masks',"Models","chromosome",'Trait'],
  right_on=['Name_backman','Masks',"Models","chromosome",'Trait']
)


pipeline_regeneron.to_csv(
  os.path.join(outdir,'pipeline_regeneron_shared.tsv.gz'),
  sep="\t",quoting=csv.QUOTE_NONE,index=False
)

pipeline_alpha.to_csv(
  os.path.join(outdir,'pipeline_alpha_shared.tsv.gz'),
  sep="\t",quoting=csv.QUOTE_NONE,index=False
)

pipeline_results.to_csv(
  os.path.join(outdir,'pipeline_results.tsv.gz'),
  sep="\t",quoting=csv.QUOTE_NONE,index=False
)

alphamissense_results.to_csv(
  os.path.join(outdir,'alphamissense_results.tsv.gz'),
  sep="\t",quoting=csv.QUOTE_NONE,index=False
)

backman_results.to_csv(
  os.path.join(outdir,'backman_results.tsv.gz'),
  sep="\t",quoting=csv.QUOTE_NONE,index=False
)



