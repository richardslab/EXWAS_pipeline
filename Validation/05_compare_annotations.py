"""
  Compare annotations
"""
import os,shutil,yaml,re,glob,csv,sqlite3
import pandas as pd


alphamiss_plof_annotations = pd.read_csv(
  "/scratch/richards/yiheng.chen/project14_ExWAS_AlphaMissense/data/Annotation/regenie.anno.file.pLOF.txt",sep="\t",header=None,quoting=csv.QUOTE_NONE
).set_axis(['snp','gene','consequence'],axis=1).query("snp.str.startswith('1:')")

pipeline_chr1_annotation = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/pipeline_results/5_1_vep_summaries_wes_qc_chr1_sitesonly.sqlite3.db"
conn = sqlite3.connect(pipeline_chr1_annotation)
cur = conn.cursor()

plof_var = cur.execute(
  f"SELECT * FROM vep_summaries WHERE SNP in ({','.join(['?']*len(alphamiss_plof_annotations))})",
  [f"chr{x}" for x in alphamiss_plof_annotations.snp.to_list()]
).fetchall()
conn.close()

plof_var_df = pd.DataFrame(
  plof_var,columns=['SNP','location','plugin','consequence','gene']
)
plof_var_df['SNP'] = plof_var_df['SNP'].apply(
  lambda val: re.sub("chr","",val)
)
shared = pd.merge(
  plof_var_df,
  alphamiss_plof_annotations,
  'inner',left_on=['SNP','gene'],right_on=['snp','gene']
)
shared_Lof = shared.query("plugin == 'IMPACT'")
set(alphamiss_plof_annotations.snp.to_list()).difference(shared_Lof.SNP.to_list())

