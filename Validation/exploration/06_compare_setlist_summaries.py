"""
  The number of variants per gene and the number of pLoF is generally 1-1 though there is slightly more from the pipeline
"""
import os,shutil,yaml,re,csv,glob,sqlite3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import tempfile
tdir = tempfile.TemporaryDirectory()
tfile = os.path.join(tdir.name,'tmp.png')

setlist_counts = pd.read_csv(
  "/scratch/richards/yiheng.chen/project14_ExWAS_AlphaMissense/data/Annotation/set_file_ALL_with_annotation_counts.csv",
  sep=",",
  quoting=csv.QUOTE_NONNUMERIC
).drop(['Unnamed: 0'],axis=1)



all_annotation_results = glob.glob(
  "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/pipeline_results/5_1_vep_summaries_wes_qc_chr*_sitesonly.sqlite3.db"
)
all_summaries = pd.DataFrame()
for each_g in setlist_counts.V1.to_list():
  test = pd.DataFrame()
  for file in all_annotation_results:
    conn = sqlite3.connect(file)
    cur = conn.cursor()
    res = cur.execute(
      "SELECT * FROM vep_summaries WHERE gene == ?",[each_g]
    ).fetchall()
    conn.close()
    test = pd.concat(
      [
        test,
        pd.DataFrame(res)
      ],axis=0
    ).reset_index(drop=True)
  if len(test) == 0:
    continue
  test.columns = ['SNP','location','plugin','consequence','gene']
  all_summaries = pd.concat(
    [all_summaries,
    pd.DataFrame(
      {
        "N_total" : len(test['SNP'].unique()),
        "N_pLoF" : test.query("plugin == 'IMPACT' & consequence == 'HIGH'").shape[0],
        "gene":each_g
      },index=[0]
    )],axis=0
  ).reset_index(drop=True)

shared_data = pd.merge(
  setlist_counts[['V1','N_total_variants',"N_pLoF_variant"]],
  all_summaries,
  'inner',left_on=['V1'],right_on=['gene']
)

fig,ax = plt.subplots(1,2)
sns.scatterplot(
  shared_data,
  x='N_total_variants',y='N_total',ax=ax[0]
)
sns.scatterplot(
  shared_data,
  x='N_pLoF',y='N_pLoF_variant',ax=ax[1]
)
ax[0].axline((0,0),slope=1)
ax[1].axline((0,0),slope=1)
fig.savefig(tfile)