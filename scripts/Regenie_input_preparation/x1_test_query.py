import shutil,os,sys,yaml,pyreadr,re,json,gzip,sqlite3
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse
from pathlib import Path
from collections import Counter

with open("/home/richards/kevin.liang2/scratch/exwas_pipeline/config/plof_or_5in5_configs/proj_config.yml",'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
CONFIG = namedtuple("params",params.keys())(**params)

conn = sqlite3.connect("/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/plof/test/5_1_vep_summaries_wes_qc_chr22_sitesonly.sqlite3.db")
cur = conn.cursor()

regeneron_masks = CONFIG.annotation_definitions['regeneron']
anno_plof = regeneron_masks['pLoF']
anno_del5in5 = {k:v for k,v in regeneron_masks['deleterious_5_of_5'].items() if k != 'var_consequence'}
for criteria,plugin_info in anno_del5in5.item():
  criteria_satisfied = Counter()
  for plugin,plugin_criteria in plugin_info.items():
    plugin_criteria_str = f" OR ".join("plugin_consequence = ?" for _ in plugin_criteria)
    plugin_var = cur.execute(
      f"""
      SELECT distinct SNP FROM vep_summaries WHERE plugin = ? AND ({plugin_criteria_str}) AND var_consequence LIKE "%missense_variant%"
      """,tuple([plugin] + plugin_criteria)
    ).fetchall()
    for x in plugin_var:
      criteria_satisfied[x[0]] +=1 
  if criteria == "all":
    valid_snps = [x for x in criteria_satisfied if criteria_satisfied[x] == len(plugin_info)]
    
    
  






