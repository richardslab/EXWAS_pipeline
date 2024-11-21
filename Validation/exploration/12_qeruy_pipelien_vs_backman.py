import os,shutil,yaml,re,csv,glob
import pandas as pd
import numpy as np

backman_res = pd.read_csv("/home/richards/kevin.liang2/scratch/exwas_pipeline/data/backman_exwas_results/GCST90083131_buildGRCh38.tsv.gz",sep="\t",quoting=csv.QUOTE_NONE)
pipeline_res = pd.DataFrame()
for i in glob.glob("/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/plof/regeneron/Regenie_S2/8_regenie_S2_OUT_wes_qc_chr*_IRNT_DBP.regenie.gz"):
  pipeline_res = pd.concat(
    [
      pipeline_res,
      pd.read_csv(i,quoting=csv.QUOTE_NONE,sep="\t",comment='#')
    ],axis=0
  ).reset_index(drop=True)

masks = {
  "M1_LoF.singleton":"M1.singleton",
  "M1_LoF.1e-05":"M1.0001",
  "M1_LoF.0.0001":"M1.001",
  "M1_LoF.0.001":"M1.01",
  "M1_LoF.0.01":"M1.1"
}