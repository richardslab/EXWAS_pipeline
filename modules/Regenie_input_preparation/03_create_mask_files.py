"""
  This creates the mask files defined in the config files
"""

import shutil,os,sys,yaml,pyreadr,re
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse
from itertools import compress

def sanity_checks():
  """Make sure we have all the inputs we need and the setup is correct so far.

  A bunch of asserts. if all good, then passes
  """

  # annotation file exists
  expected_annotation_file = os.path.join(WDIR,f'3_{VCF_NAME}_vcf_final_annotation.txt')
  assert(
    os.path.isfile(expected_annotation_file)
  ),f"Missing annotation file {os.path.isfile(expected_annotation_file)}"

  print("Check if required plugins are in VEP annotation input")
  print ("*" * 20)
  print(f"Annotation file: {expected_annotation_file}")
  extra_columns = parse_vep_headers.get_vep_plugins(expected_annotation_file)
  print(f"Vep annotation found: {extra_columns}")
  # all plugins used for mask definition exists
  for study,study_masks in CONFIG.mask_definitions.items():
    study_plugins = set()
    for mask,annotations in study_masks.items():
      for annotation,annotation_def in annotations.items():
        for plugin,plugin_criteria in annotation_def[0].items():
          study_plugins.add(plugin)

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
    mask_file = os.path.join(study_outdir,f"{VCF_NAME}_masks.txt")
    with open(mask_file,'w') as ptr:
      for mask,mask_def in study_masks.items():
        mask_def_string = ",".join(list(mask_def.keys()))
        ptr.write(f"{mask} {mask_def_string}\n")
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
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/data/wes_qc_chr3_chr_full_final.vcf.subset.sorted.vcf.gz"
    print("TEST")


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(os.path.isfile(cargs.input_vcf)),'input vcf is missing'
  assert(cargs.wdir),'output directory missing'

  print("Creating mask files")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"input VCF: {os.path.basename(cargs.input_vcf)}")
  print(f"output dir: {cargs.wdir}")
  print("="*20)


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  VCF_NAME = os.path.basename(cargs.input_vcf)
  WDIR = cargs.wdir

  sys.path.append(CONFIG.Regenie_input_prep_scripts)
  from python_helpers.vep_helpers import parse_vep_headers
  

  main()