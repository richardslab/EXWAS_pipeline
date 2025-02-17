#!/usr/bin/env python3
""" Create the annotation file for each study based on the masks defined in the config file.


"""

import shutil,os,sys,yaml,pyreadr,re,json,gzip,sqlite3
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse
from pathlib import Path
from collections import Counter

def __obtain_annotation_var(var_consequence,annotation_def,annotation_criteria,vep_summarie_file):
  plugin_criteria_vals = []
  plugin_criteria_query_str = []
  for plugin,plugin_criteria in annotation_def.items():
    if plugin in CONFIG.CONST:
      plugin_criteria_query_str += ["plugin = ? AND plugin_consequence in ({})".format(",".join(["?"]*len(plugin_criteria)))]
      plugin_criteria_vals+=[plugin]
      plugin_criteria_vals += plugin_criteria
    elif plugin in CONFIG.CONST_NUMERIC:
      threshold = float(re.sub(">|>=|<|<=","",plugin_criteria))
      if re.match("^<[0-9]+",plugin_criteria):
        plugin_criteria_query_str += ["plugin = ? AND  plugin_consequence < ?"]
        plugin_criteria_vals+=[plugin,threshold]
      elif re.match("^<=[0-9]+",plugin_criteria):
        plugin_criteria_query_str += ["plugin = ? AND  plugin_consequence <= ?"]
        plugin_criteria_vals+=[plugin,threshold]
      elif re.match("^>[0-9]+",plugin_criteria):
        plugin_criteria_query_str += ["plugin = ? AND  plugin_consequence > ?"]
        plugin_criteria_vals+=[plugin,threshold]
      elif re.match("^>=[0-9]+",plugin_criteria):
        plugin_criteria_query_str += ["plugin = ? AND  plugin_consequence >= ?"]
        plugin_criteria_vals+=[plugin,threshold]
  plugin_criteria_query_str = " OR ".join([f"({x})" for x in plugin_criteria_query_str])
  if var_consequence is not None:
    var_conseq_str = " OR ".join("var_consequence LIKE ?" for _ in var_consequence)
    var_consequence = [f"%{x}%" for x in var_consequence]
    try:
      conn = sqlite3.connect(vep_summarie_file)
      cur = conn.cursor()
      plugin_var = cur.execute(
        f"""SELECT SNP,gene FROM vep_summaries WHERE ({plugin_criteria_query_str}) AND ({var_conseq_str})""",tuple(plugin_criteria_vals + var_consequence)
      ).fetchall()
      conn.close()
    except Exception as e:
      print(f"SQLITE3 error fetching pairs: {e}")
      conn.close()
      raise
  else:
    try:
      conn = sqlite3.connect(vep_summarie_file)
      cur = conn.cursor()
      plugin_var = cur.execute(
        f"""SELECT SNP,gene FROM vep_summaries WHERE plugin = ? AND ({plugin_criteria_query_str})
        """,tuple([plugin] + plugin_criteria_vals)
      ).fetchall()
      conn.close()
    except Exception as e:
      print(f"SQLITE3 error fetching pairs: {e}")
      conn.close()
      raise
  criteria_counter = Counter()
  for x in plugin_var:
    criteria_counter[x] += 1
  satisfied_var = []
  if annotation_criteria == "all":
    satisfied_var = [x for x in criteria_counter if criteria_counter[x] == len(annotation_def)]
  elif annotation_criteria == "any":
    satisfied_var = list(criteria_counter.values())
  else:
    satisfied_var = [x for x in criteria_counter if criteria_counter[x] >= annotation_criteria]
  return satisfied_var


def main():
  print("creating annotation file per study")
  print("*"*20)
  vep_summarie_file = os.path.join(
    ANNOTATION_DIR,f"5_1_vep_summaries_{VCF_NAME}.sqlite3.db"
  )
  assert(
    os.path.isfile(vep_summarie_file)
  ),f"missing summary file"
  
  # Make sure there isn't an annotation file already
  all_studies = list(CONFIG.mask_definitions.keys())
  for study in all_studies:
    os.makedirs(os.path.join(WDIR,study),exist_ok=True)
    study_ofile = os.path.join(WDIR,study,f"annotations_{VCF_NAME}.txt")
    assert(
      not os.path.isfile(study_ofile)
    ),f"annotation file found for {study}. Delete it first, or skip this step"
    study_annotation_info = CONFIG.annotation_definitions[study]
    study_ofile = os.path.join(WDIR.rstrip("/"),study,f"annotations_{VCF_NAME}.txt")
    study_var = dict()
    with open(study_ofile,'w') as ptr:
      for annotation in CONFIG.annotation_order[study]:
        annotation_info = CONFIG.annotation_definitions[study][annotation]
        all_annotations = {k:v for k,v in annotation_info.items() if k != 'var_consequence'}
        if 'var_consequence' in annotation_info:
          var_consequence = annotation_info['var_consequence']
        else:
          var_consequence = None
        all_var = []
        for annotation_criteria,annotation_def in all_annotations.items():
          annotation_var = __obtain_annotation_var(var_consequence,annotation_def,annotation_criteria,vep_summarie_file)
          all_var += annotation_var
          all_var = list(set(all_var))
        for var in all_var:
          if var not in study_var:
            study_var[var] = annotation
            res = ptr.write(f"{var[0]}\t{var[1]}\t{annotation}\n")
  return

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--config_file','-c',
    dest='cfile',
    help='configuration yaml file'
  )
  parser.add_argument(
    '--input_vcf','-i',
    dest='input_vcf',
    help="input VCF file",
    type=str
  )
  parser.add_argument(
    '--annotation_summary_dir',
    dest='anno_dir',
    help="input directory",
    type=str
  )
  parser.add_argument(
    '--test',
    default='f',
    type=str
  )
  cargs =   parser.parse_args()

  if cargs.test =='t':
    from unittest import mock
    cargs = mock.Mock()
    cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/alphamis_exact_configs/proj_config.yml"
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/alphamissense_results/ukb_merged_1-22_sitesonly.vcf"
    __file__ = "/home/richards/kevin.liang2/scratch/exwas_pipeline/src/modules/Regenie_input_preparation/04_2_create_annotation_files.py"
    print("TEST")

  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)

  VCF_NAME = Path(cargs.input_vcf).stem
  ANNOTATION_DIR = cargs.anno_dir
  WDIR=os.getcwd()

  print("Creating annotation file")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"input VCF: {os.path.basename(cargs.input_vcf)}")
  print("="*20)


  main()