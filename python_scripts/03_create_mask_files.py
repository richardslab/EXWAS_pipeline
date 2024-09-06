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
  expected_annotation_file = os.path.join(CONFIG.wdir,f'2_{VCF_NAME}_vcf_final_annotation.txt')
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

  
  # all masks have definitions
  for study,mask_names in CONFIG.mask_names.items():
    all_definitions_provided = list(CONFIG.mask_definitions[study].keys())
    study_mask_definitions_needed = []
    for x in list(mask_names.values()):
      x = x.split(",")
      study_mask_definitions_needed += [x.strip() for x in x]
    study_mask_definitions_needed = list(set(study_mask_definitions_needed))
    assert(
      all(
        [x in all_definitions_provided for x in study_mask_definitions_needed]
      )
    ),f"{study}: not all definitions are provided for masks"
      
  
  return


def main():
  sanity_checks()
  # write the mask file for each study
  for study,mask_names in CONFIG.mask_names.items():

    study_outdir = os.path.join(CONFIG.wdir,study)
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
  cargs =   parser.parse_args()

  import mock
  cargs = mock.Mock()
  cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  print(f"Using {os.path.basename(cargs.cfile)}")


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)

  sys.path.append(CONFIG.script_dir)
  from python_scripts.python_helpers.vep_helpers import parse_vep_headers
  
  VCF_NAME = os.path.basename(CONFIG.input_vcf)

  main()