""" Create the annotation file for each study based on the masks defined in the config file.


"""

import shutil,os,sys,yaml,pyreadr,re,json,gzip,sqlite3
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse
from pathlib import Path

def __init_annotation():
  all_studies = list(CONFIG.mask_definitions.keys())
  for study in all_studies:
    study_ofile = os.path.join(WDIR,study,f"annotations_{VCF_NAME}.txt")
    assert(
      not os.path.isfile(study_ofile)
    ),f"annotation file found for {study}. Delete it first, or skip this step"    
    study_ofile = Path(study_ofile)
    study_ofile.touch(exist_ok=True)
    print(f"annotation for {study} written to {study_ofile}")
  return all_studies

def __eval_plugins(var,gene,vep_summarie_file,criteria,plugin_definition):
  n_criteria_satisfied = 0
  for plugin,plugin_criteria in plugin_definition.items():
    try:
      conn2 = sqlite3.connect(vep_summarie_file)
      cur2 = conn2.cursor()
      variant_criteria = cur2.execute(
        """
        SELECT * FROM vep_summaries WHERE SNP = :SNP AND gene = :gene AND plugin = :plugin
        """,{"SNP":var,"gene":gene,"plugin":plugin}
      ).fetchall()
      conn2.close()
    except Exception as e:
      print(f"SQLITE3 error fetching annoations: {e}")
      conn2.close()
      raise
    if len(variant_criteria) == 0:
      continue
    assert(len(variant_criteria)==1),f"{var} {gene} {plugin} had more than 1 entry for gene plugin"
    variant_criteria = variant_criteria[0]
    if plugin in CONFIG.CONST:
      if variant_criteria[3] in plugin_criteria:
        n_criteria_satisfied += 1
        continue
    elif plugin in CONFIG.CONST_NUMERIC:
      if re.match(">|>=",plugin_criteria):
        threshold = float(re.sub(">|>=","",plugin_criteria))
        if re.match(">",plugin_criteria) and variant_criteria[3] > threshold:
          n_criteria_satisfied += 1
          continue
        elif re.match(">=",plugin_criteria) and variant_criteria[3] >= threshold:
          n_criteria_satisfied+=1
          continue
        else:
          assert(False),f"Parse plugin error: plugin {plugin}, criteria {plugin_criteria}"
      elif re.match("<|<=",plugin_criteria):
        threshold = float(re.sub("<|<=",plugin_criteria))
        if re.match("<",plugin_criteria) and variant_criteria[3] < threshold:
          n_criteria_satisfied+=1
          continue
        elif re.match("<=",plugin_criteria) and variant_criteria[3] <= threshold:
          n_criteria_satisfied+=1
          continue
        else:
          assert(False),f"Parse plugin error: plugin {plugin}, criteria {plugin_criteria}"
      else:
        assert(False),f"Parse plugin error: plugin {plugin}, criteria {plugin_criteria}"
  if criteria == "all":
    if n_criteria_satisfied == len(plugin_definition):
      return True
    else:
      return False
  elif criteria == "any":
    if n_criteria_satisfied >= 1:
      return True
    else:
      return False
  else:
    criteria = int(criteria)
    if n_criteria_satisfied >= criteria:
      return True
    else:
      return False


def __eval_annotation(var,gene,annotation_def,vep_summarie_file):
  satisfied = []
  for criteria,plugin_definition in annotation_def.items():
    plugin_status = __eval_plugins(var,gene,vep_summarie_file,criteria,plugin_definition)
    satisfied.append(plugin_status)
  if all(satisfied):
    return True
  else:
    return False
  
  

def main():
  print("creating annotation file per study")
  print("*"*20)
  vep_summarie_file = os.path.join(
    WDIR,f"5_1_vep_summaries_{VCF_NAME}.sqlite3.db"
  )
  assert(
    os.path.isfile(vep_summarie_file)
  ),f"missing summary file"

  try:
    conn = sqlite3.connect(vep_summarie_file)
    cur = conn.cursor()
    n_rows = cur.execute(
      """
      SELECT count(*) FROM (SELECT distinct SNP,gene FROM vep_summaries)
      """
    ).fetchone()
    n_rows = float(n_rows[0])
    variant_genes = cur.execute(
      """
      SELECT distinct SNP,gene FROM vep_summaries
      """
    )
  except Exception as e:
    print(f"SQLITE3 error fetching pairs: {e}")
    conn.close()
    raise
  
  # Make sure there isn't an annotation file created already
  all_studies = __init_annotation()

  # it loops through all variants for each gene
  line_ct = 0
  for var,gene in variant_genes:
    line_ct += 1
    if line_ct % 5000 == 0:
      print(f"{line_ct}/{n_rows}: {'{:.3e}'.format(line_ct/n_rows*100)}% completed")
    for study in all_studies:
      study_annotations = CONFIG.annotation_definitions[study]
      # assign this variant to masks
      var_with_possible_annotations = []
      for annotation,annotation_def in study_annotations.items():
        anno_fit = __eval_annotation(var,gene,annotation_def,vep_summarie_file)
        if anno_fit:
          var_with_possible_annotations.append(annotation)
          
      # assign variant to a mask based on order
      annotation_order = CONFIG.annotation_order[study]
      if len(var_with_possible_annotations) > 0:
        var_anno = None
        for i in annotation_order:
          if i in var_with_possible_annotations:
            var_anno = i
            break
        assert(var_anno != None),f"did not assign unique mask {gene} {var} {study} {var_with_possible_annotations}"
        study_ofile = os.path.join(WDIR,study,f"annotations_{VCF_NAME}.txt")
        with open(study_ofile,'a') as ptr:
          res = ptr.write(f"{var}\t{gene}\t{var_anno}\n")
  conn.close()
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
    '--wdir',
    dest='wdir',
    help="Output directory",
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
    cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"
    cargs.wdir="/scratch/richards/kevin.liang2/exwas_pipeline/results/pipeline_results"
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/sitesonly_VCF/wes_qc_chr2_sitesonly.vcf"
    __file__ = "/home/richards/kevin.liang2/scratch/exwas_pipeline/src/modules/Regenie_input_preparation/04_1_create_annotation_summaries.py"
    print("TEST")

  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)

  VCF_NAME = Path(cargs.input_vcf).stem
  WDIR = cargs.wdir

  print("Creating annotation file")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"input VCF: {os.path.basename(cargs.input_vcf)}")
  print(f"output dir: {cargs.wdir}")
  print("="*20)


  main()