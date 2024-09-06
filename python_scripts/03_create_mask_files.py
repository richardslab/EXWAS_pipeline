"""
  This creates the mask files defined in the config files
"""

import shutil,os,sys,yaml,pyreadr,re
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse

def sanity_checks():
  """Make sure we have all the inputs we need and the setup is correct so far
  """

  # annotation file
  expected_annotation_file = os.path.join(CONFIG.wdir,'2_vcf_final_annotation.txt')
  assert(
    os.path.isfile(expected_annotation_file)
  ),f"Missing annotation file {os.path.isfile(expected_annotation_file)}"

  # the maks definitions are found within the annotation file
  ## obtain all the plugins that was attempted by extracting all 'extra' columns in the vep outputs

  # the extra column key is the last set of info before the results
  # https://useast.ensembl.org/info/docs/tools/vep/vep_formats.html
  start_serach_str = "## Extra column keys:"
  end_search_str="#Uploaded_variation"

  extra_columns = []
  with open(expected_annotation_file,'r') as ptr:
    start_search = False
    end_search = False
    for line in ptr:
      if re.match(start_serach_str,line):
        start_search=True
        continue
      if re.match(end_search_str,line):
        end_search = True
        break
      if start_search and not end_search:
        extra_columns += [re.sub("##","",line.split(":")[0]).strip()]
  
  # Make sure the plugins used for mask definitions
  # are found
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

  
  # make sure we have definition for all the mask filter
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
      
  # if all assertion passes, it's all good
  return


def main():
  sanity_checks()


  return


if __name__ == "__main__":
  # parser = argparse.ArgumentParser()
  # parser.add_argument(
  #   '--config_file','-c',
  #   dest='cfile',
  #   default="/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml",
  #   help='configuration yaml file'
  # )
  # cargs =   parser.parse_args()

  import mock
  cargs = mock.Mock()
  cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  print(f"Using {os.path.basename(cargs.cfile)}")


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)

  main()