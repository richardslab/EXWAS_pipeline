""" Create the annotation file for each study based on the masks defined in the config file.

"""

import shutil,os,sys,yaml,pyreadr,re,json,gzip
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse

def main():
  vep_summarie_file = os.path.join(
    CONFIG.wdir,'vep_consequence_summaries',f"{VCF_NAME}_vep_summaries.json.gz"
  )
  assert(
    os.path.isfile(vep_summarie_file)
  ),f"missing summary file"

  with gzip.open(vep_summarie_file,'r') as ptr:
    vep_summaries = json.load(ptr)

  all_studies = list(CONFIG.mask_names.keys())

  for gene,var in vep_summaries.items():
    for study in all_studies:
      study_annotations = CONFIG.mask_definition[study]


  


  return 

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--config_file','-c',
    dest='cfile',
    default="/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml",
    help='configuration yaml file'
  )
  cargs =   parser.parse_args()

  # import mock
  # cargs = mock.Mock()
  # cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  print(f"Using {cargs.cfile}")



  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)


  CONSTANT = CONFIG.CONST
  CONST_NUMERIC = CONFIG.CONST_NUMERIC
  

  VCF_NAME = os.path.basename(CONFIG.input_vcf)
  sys.path.append(CONFIG.script_dir)




  main()