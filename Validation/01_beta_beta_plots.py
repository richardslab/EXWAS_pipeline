"""
  Compare the effect sizes my runs and Backman et al 2021
"""
#%%

import os,shutil,yaml,pickle,re,glob
import pandas as pd
import numpy as np
import matplolib.pyplot as plt
import seaborn as sns

#%%
# downloaded_data_constants
backman_gwas_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/data/backman_exwas_results"
backman_gwas_files = {
  "SBP":"GCST90083132_buildGRCh38.tsv.gz",
  "DBP":"GCST90083131_buildGRCh38.tsv.gz",
  "Standing_height":"GCST90083263_buildGRCh38.tsv.gz",
  "LDL":"GCST90083019_buildGRCh38.tsv.gz",
  "TG":"GCST90083030_buildGRCh38.tsv.gz",
  "Calcium":"GCST90083009_buildGRCh38.tsv.gz",
  "Dbilirubin":"GCST90083007_buildGRCh38.tsv.gz",
  "Glucose":"GCST90083015_buildGRCh38.tsv.gz",
  "RBC":"GCST90082965_buildGRCh38.tsv.gz",
  "BMI":"GCST90082670_buildGRCh38.tsv.gz"
}
# own regenie results:
pipeline_result_path="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/PAST_pipeline_results/intial_runs/regeneron/Regenie_S2"
pipeline_result_files = {
  "SBP":"8_regenie_S2_OUT_wes_qc_chr10_SBP.regenie.gz",
  "DBP":"8_regenie_S2_OUT_wes_qc_chr10_DBP.regenie.gz",
  "Standing_height":"8_regenie_S2_OUT_wes_qc_chr10_Standing_height.regenie.gz",
  "LDL":"8_regenie_S2_OUT_wes_qc_chr10_LDL.regenie.gz",
  "TG":"8_regenie_S2_OUT_wes_qc_chr10_TG.regenie.gz",
  "Calcium":"8_regenie_S2_OUT_wes_qc_chr10_Calcium.regenie.gz",
  "Dbilirubin":"8_regenie_S2_OUT_wes_qc_chr10_Dbilirubin.regenie.gz",
  "Glucose":"8_regenie_S2_OUT_wes_qc_chr10_Glucose.regenie.gz",
  "RBC":"8_regenie_S2_OUT_wes_qc_chr10_RBC.regenie.gz",
  "BMI":"8_regenie_S2_OUT_wes_qc_chr10_BMI.regenie.gz",
  "HDL":"8_regenie_S2_OUT_wes_qc_chr10_HDL.regenie.gz"
}

# Compare results

