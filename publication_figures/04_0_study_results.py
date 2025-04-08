"""
Compare significant results between studies

gather study results
"""

import os,shutil,yaml,re,csv,scipy,glob
import pandas as pd
import numpy as np

with open("./parameters/script_params.yml",'r') as ptr:
  params = yaml.safe_load(ptr)['study_results_0']

outdir = params['outdir']

#%% all backman results
backman_gwas_path=params['backman_gwas_path']
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

backman_masks = {
  "M1.singleton" : "M1_LoF.singleton",
  "M1.0001" : "M1_LoF.1e-05",
  "M1.001" : "M1_LoF.0.0001",
  "M1.01" : "M1_LoF.0.001",
  "M1.1" : "M1_LoF.0.01",
  'M3.singleton' : "M3_LoF_or_deleterious_5_of_5.singleton",
  'M3.0001' : "M3_LoF_or_deleterious_5_of_5.1e-05",
  'M3.001' : "M3_LoF_or_deleterious_5_of_5.0.0001",
  'M3.01' : "M3_LoF_or_deleterious_5_of_5.0.001",
  'M3.1' : "M3_LoF_or_deleterious_5_of_5.0.01"
}
all_backman_results = pd.DataFrame()
for pheno,file in backman_gwas_files.items():
  backman_res = pd.read_csv(
    os.path.join(backman_gwas_path,file),sep="\t",quoting=csv.QUOTE_NONE
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
  )
  backman_res['Trait'] = pheno
  all_backman_results = pd.concat(
    [all_backman_results,backman_res],
    axis=0
  ).reset_index(drop=True)


#%% alphamissense results
alphamiss_gwas_path = params['alphamiss_gwas_path']
alphamiss_continous_gwas_files = {"ZBMD":"step_2_ZBMD_ExWAS_pLOF_missense_ZBMD.regenie",
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
alphamiss_binary_gwas_files = {
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

Alphamissense_masks = {
  "pLOF_only.singleton": "M1_LoF.singleton",
  "pLOF_only.0.01": "M1_LoF.0.01",
  "pLOF_only.0.001": "M1_LoF.0.001",
  "pLOF_and_55missense.singleton": 'M2_LoF_or_deleterious.singleton',
  "pLOF_and_55missense.0.01": 'M2_LoF_or_deleterious.0.01',
  "pLOF_and_55missense.0.001": 'M2_LoF_or_deleterious.0.001'
}

test_continuous_map = {
  "ADD":'ADD-WGR-LR', 
  "ADD-SKAT":'ADD-WGR-SKAT', 
  "ADD-SKATO":'ADD-WGR-SKATO',
  "ADD-BURDEN-ACAT":'ADD-WGR-BURDEN-ACAT'
}


all_alphamiss_res = pd.DataFrame()
for pheno,file in alphamiss_continous_gwas_files.items():
  Alphamissense_res = pd.read_csv(
    os.path.join(alphamiss_gwas_path,file),
    sep=" ",quoting=csv.QUOTE_NONE,comment="#"
  ).assign(
    Masks = lambda df: df['ALLELE1'].map(Alphamissense_masks),
    p_value = lambda df: df['LOG10P'].apply(lambda val: 10**(-1 * val)),
    Models = lambda df: df['TEST'].map(test_continuous_map)
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
  )
  Alphamissense_res['Trait'] = pheno
  all_alphamiss_res = pd.concat(
    [all_alphamiss_res,Alphamissense_res],
    axis=0
  ).reset_index(drop=True)

test_binary_map = {
  "ADD":'ADD-WGR-FIRTH', 
  "ADD-SKAT":'ADD-WGR-SKAT', 
  "ADD-SKATO":'ADD-WGR-SKATO',
  "ADD-BURDEN-ACAT":'ADD-WGR-BURDEN-ACAT'
}
for pheno,file in alphamiss_binary_gwas_files.items():
  Alphamissense_res = pd.read_csv(
    os.path.join(alphamiss_gwas_path,file),
    sep=" ",quoting=csv.QUOTE_NONE,comment="#"
  ).assign(
    Masks = lambda df: df['ALLELE1'].map(Alphamissense_masks),
    p_value = lambda df: df['LOG10P'].apply(lambda val: 10**(-1 * val)),
    Models = lambda df: df['TEST'].map(test_binary_map)
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
  )
  Alphamissense_res['Trait'] = pheno
  all_alphamiss_res = pd.concat(
    [all_alphamiss_res,Alphamissense_res],
    axis=0
  ).reset_index(drop=True)

# save
all_backman_results.to_csv(
  os.path.join(outdir,"all_backman_results.tsv.gz"),
  sep="\t",quoting=csv.QUOTE_NONE,index=False
)
all_alphamiss_res.to_csv(
  os.path.join(outdir,"all_alphamiss_res.tsv.gz"),
  sep="\t",quoting=csv.QUOTE_NONE,index=False
)