#!/usr/bin/env python3
"""
  This creates the mask files defined in the config files
"""

import shutil,os,sys,yaml,pyreadr,re
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse
from itertools import compress
from pathlib import Path

def sanity_checks():
  """Make sure we have all the inputs we need and the setup is correct so far.

  A bunch of asserts. if all good, then passes
  """
  # annotation file exists
  expected_annotation_file = os.path.join(IDIR,f'3_annotation_results_{VCF_NAME}.txt.gz')
  assert(
    os.path.isfile(expected_annotation_file)
  ),f"Missing annotation file {expected_annotation_file}"
  print("Check if required plugins are in VEP annotation input")
  print ("*" * 20)
  print(f"Annotation file: {expected_annotation_file}")
  extra_columns = parse_vep_headers.get_vep_plugins(expected_annotation_file)
  print(f"Vep annotation found: {extra_columns}")
  # all plugins used for mask definition exists
  for study,annotation_information in CONFIG.annotation_definitions.items():
    study_plugins = set()
    for annotation_name,plugin_information in annotation_information.items():
      plugin_information = {k:v for k,v in plugin_information.items() if k!= 'var_consequence'}
      for criteria,plugin_criteria in plugin_information.items():
        study_plugins = study_plugins.union(set(list(plugin_criteria.keys())))
    plugin_membership = [x in extra_columns for x in study_plugins]
    plugin_missing_val = list(compress(
      study_plugins,[not x for x in plugin_membership]
    ))
    assert(
      all(
        plugin_membership
      )
    ),f"Plugins needed for {study} not found: {plugin_missing_val}"
  print("Required plugins are in annotation file")
  print("="*20)
  return


def main():
  sanity_checks()

  # write the mask file for each study
  for study,study_masks in CONFIG.mask_definitions.items():
    study_outdir = os.path.join(WDIR,study)
    os.makedirs(study_outdir,exist_ok=True)
    mask_file = os.path.join(study_outdir,f"masks_{VCF_NAME}.txt")
    with open(mask_file,'w') as ptr:
      for mask_name,mask_def in study_masks.items():
        mask_def_string = ",".join(mask_def)
        res = ptr.write(f"{mask_name} {mask_def_string}\n")
    print(f"Masks for {study} written to {study_outdir}")
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
    '--input_dir',
    dest='idir',
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
    cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/sitesonly_VCF/wes_qc_chr10_sitesonly.vcf"
    __file__ = "/home/richards/kevin.liang2/scratch/exwas_pipeline/src/modules/Regenie_input_preparation/03_create_mask_files.py"
    print("TEST")


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(os.path.isfile(cargs.input_vcf)),'input vcf is missing'

  print("Creating mask files")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"input VCF: {os.path.basename(cargs.input_vcf)}")
  print("="*20)


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  VCF_NAME = Path(cargs.input_vcf).stem
  WDIR = os.getcwd()
  IDIR = cargs.idir

  sys.path.append(os.path.dirname(__file__))
  from python_helpers.vep_helpers import parse_vep_headers
  

  main()