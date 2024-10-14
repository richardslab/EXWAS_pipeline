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
    study_ofile = os.path.join(WDIR,study,f"{VCF_NAME}_annotations.txt")
    assert(
      not os.path.isfile(study_ofile)
    ),f"annotation file found for {study}. Delete it first, or skip this step"

    study_ofile = Path(study_ofile)
    study_ofile.touch(exist_ok=True)
    print(f"annotation for {study} written to {study_ofile}")
  return all_studies

def __obtain_study_annotations(study):
  study_masks = CONFIG.mask_definitions[study]
  all_annotation_info = dict()
  for mask_def in study_masks.values():
    all_annotation_info.update(mask_def)
  return all_annotation_info

def __eval_annotation(var,gene,annotation_def,cur):
  annotation_type = annotation_def[1]
  annotation_criteria = annotation_def[0]
  n_criteria_satisfied = 0
  for plugin,plugin_criteria in annotation_criteria.items():
    variant_criteria = cur.execute(
      """
      SELECT * FROM vep_summaries WHERE SNP = :SNP AND gene = :gene AND plugin = :plugin
      """,{"SNP":var,"gene":gene,"plugin":plugin}
    ).fetchall()
    if len(variant_criteria) == 0:
      continue
    assert(len(variant_criteria)==1),f"{var} {gene} {plugin} had more than 1 entry"
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
  
  if annotation_type == "all":
    if n_criteria_satisfied == len(annotation_criteria):
      return True
    else:
      return False
  elif annotation_type == "any":
    if n_criteria_satisfied >= 1:
      return True
    else:
      return False
  else:
    annotation_type = int(annotation_type)
    if n_criteria_satisfied >= annotation_type:
      return True
    else:
      return False
  

def main():
  print("creating annotation file per study")
  print("*"*20)
  vep_summarie_file = os.path.join(
    WDIR,f"5_1_{VCF_NAME}_vep_summaries.sqlite3.db"
  )
  assert(
    os.path.isfile(vep_summarie_file)
  ),f"missing summary file"

  conn = sqlite3.connect(vep_summarie_file)
  cur = conn.cursor()
  variant_genes = cur.execute(
    """
    SELECT distinct SNP,gene FROM vep_summaries
    """
  )


  # Make sure there isn't an annotation file created already
  all_studies = __init_annotation()


  # it loops through all variants for each gene
  for var,gene in variant_genes:
    for study in all_studies:
      study_annotations = __obtain_study_annotations(study)
      # assign this variant to masks
      var_with_possible_annotations = []
      for annotation,annotation_def in study_annotations.items():
        anno_fit = __eval_annotation(var,gene,annotation_def,cur)
        if anno_fit:
          var_with_possible_annotations.append(annotation)
          
      # assign variant to a mask based on order
      annotation_order = CONFIG.annotation_order[study]
      if len(var_with_possible_annotations) > 0:
        var_anno = None
        for i in var_with_possible_annotations:
          if i in annotation_order:
            var_anno = i
            break
        assert(var_anno != None),f"did not assign unique mask {gene} {var} {study} {var_with_possible_annotations}"
        study_ofile = os.path.join(WDIR,study,f"{VCF_NAME}_annotations.txt")
        with open(study_ofile,'a') as ptr:
          ptr.write(
            f"{var}\t{gene}\t{var_anno}\n"
          )
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
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/data/example_homo_sapiens_GRCh38.vcf"
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