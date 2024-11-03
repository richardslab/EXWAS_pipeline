# Majority (98%) of pLOF annotated by pipeline is in alphamis annotated as pLoF
# Many (25%) annotated as plof alphamis not in pipeline
# # these variants seems was dropped by vep?
# # if true, then 200K variants are dropped...

import sys,os,re,glob,shutil,yaml,csv
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


alphamiss_plof = pd.read_csv(
  "/scratch/richards/yiheng.chen/project14_ExWAS_AlphaMissense/data/Annotation/regenie.anno.file.pLOF.txt",
  sep="\t",quoting=csv.QUOTE_NONE,header=None
).set_axis(['snp','gene','consequence'],axis=1).query("snp.str.startswith('2:')")

pipeline_plof = pd.DataFrame()
for f in glob.glob("/home/richards/kevin.liang2/scratch/exwas_pipeline/results/pipeline_results/regeneron/annotations_wes_qc_chr*_sitesonly.txt"):
  pipeline_chr = pd.read_csv(
    f,
    sep="\t",header=None,quoting=csv.QUOTE_NONE
  ).set_axis(['snp','gene','consequence'],axis=1)
  pipeline_chr['snp'] = pipeline_chr['snp'].apply(
    lambda val: re.sub("chr","",val)
  )
  pipeline_plof = pd.concat([pipeline_plof,pipeline_chr],axis=0).reset_index(drop=True)

unique_pipeline_plof = pipeline_plof.snp.drop_duplicates().dropna().to_list()
unique_alphamis_plof = alphamiss_plof.snp.drop_duplicates().dropna().to_list()
shared_snps = set(unique_alphamis_plof).intersection(
  set(unique_pipeline_plof)
)
len(shared_snps)/len(unique_pipeline_plof)
len(shared_snps)/len(unique_alphamis_plof)


alphamis_only = set(unique_alphamis_plof).difference(unique_pipeline_plof)
len(alphamis_only)
