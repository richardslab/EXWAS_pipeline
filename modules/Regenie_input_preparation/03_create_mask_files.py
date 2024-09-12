"""
  This creates the mask files defined in the config files
"""

import shutil,os,sys,yaml,pyreadr,re
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse

def sanity_checks():
  """Make sure we have all the inputs we need and the setup is correct so far.

  A bunch of asserts. if all good, then passes
  """

  # annotation file exists
  expected_annotation_file = os.path.join(WDIR,f'2_{VCF_NAME}_vcf_final_annotation.txt')
  assert(
    os.path.isfile(expected_annotation_file)
  ),f"Missing annotation file {os.path.isfile(expected_annotation_file)}"

  
  extra_columns = parse_vep_headers.get_vep_plugins(expected_annotation_file)
  # all plugins used for mask definition exists
  for study in CONFIG.mask_definitions.keys():
    study_masks = CONFIG.mask_definitions[study]
    study_plugins = []
    for mask_def in study_masks.values():
      study_plugins += list(set(mask_def.keys()))
    assert(
      all(
        [x in extra_columns for x in study_plugins]
      )
    ),f"Plugins needed for {study} not found"

  
  
  return


def main():
  sanity_checks()
  # write the mask file for each study
  for study,mask_names in CONFIG.mask_names.items():

    study_outdir = os.path.join(WDIR,study)
    os.makedirs(study_outdir,exist_ok=True)

    mask_file = os.path.join(study_outdir,f"{study}_masks.txt")
    with open(mask_file,'w') as ptr:
      for mask_name,mask_def in mask_names.items():
        ptr.write(f"{mask_name} {mask_def}\n")


  return


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--config_file','-c',
    dest='cfile',
    default="/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml",
    help='configuration yaml file'
  )
  parser.add_argument(
    '--input_vcf','-i',
    dest='input_vcf',
    nargs=1,
    help="input VCF file",
    type=str
  )
  parser.add_argument(
    '--wdir',
    dest='wdir',
    nargs=1,
    help="Output directory",
    type=str
  )
  cargs =   parser.parse_args()


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(os.path.isfile(cargs.input_vcf)),'input vcf is missing'
  assert(cargs.wdir),'output directory missing'
  print(f"Using {os.path.basename(cargs.cfile)}")
  print(f"Using {os.path.basename(cargs.input_vcf)}")
  print(f"Outputs in {cargs.wdir}")


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  VCF_NAME = os.path.basename(cargs.input_vcf)
  WDIR = cargs.wdir

  sys.path.append(CONFIG.script_dir)
  from python_scripts.python_helpers.vep_helpers import parse_vep_headers
  

  main()