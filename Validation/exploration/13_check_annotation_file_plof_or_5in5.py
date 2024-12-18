import os,shutil,yaml,re,csv,scipy,glob,sqlite3
import pandas as pd
import numpy as np

pipeline_plof_or_del5in5 = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/alphamiss_exact/alphamiss_plof_5in5/annotations_ukb_merged_1-22_sitesonly.txt"

alpha_plof_or_del5in5 = pd.read_csv(
  "/scratch/richards/yiheng.chen/project14_ExWAS_AlphaMissense/data/Annotation/annotation_files_used_for_ExWAS/regenie.anno.file.pLOF.missense.txt",sep="\t",quoting=csv.QUOTE_NONE,header=None
)

alpha_plof_or_del5in5.columns = [f"alpha_{x}" for x in ["SNP","GENE","annotation"]]
annotation_map = {
  "pLoF" : "pLoF",
  "missense.1in5" : "missense.1in5",
  "missense.5in5" : 'deleterious_5_of_5'
}

alpha_plof_or_del5in5['alpha_annotation_matched'] = alpha_plof_or_del5in5["alpha_annotation"].map(annotation_map)

alpha_plof_or_del5in5 = alpha_plof_or_del5in5.assign(
  iid = lambda df: df['alpha_SNP'] + "-"+df['alpha_annotation_matched']
)

pipeline_plof_or_del5in5_annotation = pd.read_csv(
  pipeline_plof_or_del5in5,sep="\t",header=None,quoting=csv.QUOTE_NONE
)
pipeline_plof_or_del5in5_annotation.columns = [f"pipeline_{x}" for x in ["SNP","GENE","annotation"]]

pipeline_plof_or_del5in5_annotation = pipeline_plof_or_del5in5_annotation.assign(
  id = lambda df: df['pipeline_SNP']+"-" + df["pipeline_annotation"]
)



# in alpha not in pipeline
alpha_only_vars = set(alpha_plof_or_del5in5.alpha_SNP).difference(
  set(pipeline_plof_or_del5in5_annotation.pipeline_SNP)
)


alpha_only_var_annotatinos = alpha_plof_or_del5in5.query(
  "alpha_SNP.isin(@alpha_only_vars)"
) 
alpha_only_var_annotatinos.alpha_annotation.unique()

alpha_only_var_with_plof = alpha_only_var_annotatinos.query("alpha_annotation == 'pLoF'") 

annotation_summaries = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/alphamiss_exact/testing/testing_annotation_sqlite3_updated.db"

# vep annotation for these variants
conn = sqlite3.connect(annotation_summaries)
cur = conn.cursor()
alpha_only_plof_var_annotations = pd.DataFrame()
for var in alpha_only_var_with_plof.alpha_SNP.to_list():  
  query = "SELECT SNP,plugin,plugin_consequence FROM vep_summaries WHERE SNP = ? AND plugin = 'IMPACT'"
  res = cur.execute(
    query,[var]).fetchall()
  if len(res) == 0:
    alpha_only_plof_var_annotations = pd.concat(
      [
        alpha_only_plof_var_annotations,
        pd.DataFrame(
          {
            "var":var,
            "plugin":None,
            "plugin_consequence":None
          },index=[0]
        )
      ],axis=0
    ).reset_index(drop=True)  
  else:
    alpha_only_plof_var_annotations = pd.concat(
      [
        alpha_only_plof_var_annotations,
        pd.DataFrame(res).set_axis(
          ['var','plugin','plugin_consequence'],axis=1
        )
      ],axis=0
    ).reset_index(drop=True)  

conn.close()

alpha_only_plof_var_annotations.groupby(['plugin_consequence']).count()
missing_plof = alpha_only_plof_var_annotations[[x is None for x in alpha_only_plof_var_annotations.plugin]]
