""" Create the annotation file for each study based on the masks defined in the config file.

"""

import shutil,os,sys,yaml,pyreadr,re,json,gzip
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse
from pathlib import Path

def main():
  vep_summarie_file = os.path.join(
    WDIR,'vep_consequence_summaries',f"{VCF_NAME}_vep_summaries.json.gz"
  )
  assert(
    os.path.isfile(vep_summarie_file)
  ),f"missing summary file"

  with gzip.open(vep_summarie_file,'r') as ptr:
    vep_summaries = json.load(ptr)

  # Make sure there isn't an annotation file created already
  all_studies = list(CONFIG.mask_names.keys())
  for study in all_studies:
    study_ofile = os.path.join(WDIR,study,f"{VCF_NAME}_annotations.txt")
    assert(
      not os.path.isfile(study_ofile)
    ),f"annotation file found for {study}. Delete it first, or skip this step"
    study_ofile = Path(study_ofile)
    study_ofile.touch(exist_ok=True)


  # it loops through all variants for each gene
  for gene,var in vep_summaries.items():
    for var,var_info in var.items():
      for study in all_studies:
        # assign this variant to masks
        study_definitions = CONFIG.mask_definitions[study]
        var_with_possible_annotations = []
        for mask_name,plugin_defs in study_definitions.items():
          for plugin,plugin_criteria in plugin_defs.items():
            if plugin in var_info['consequences']:
              if var_info['consequences'][plugin] != None:
                if plugin in CONFIG.CONST:
                  if var_info['consequences'][plugin] in plugin_criteria:
                    var_with_possible_annotations.append(mask_name)
                elif plugin in CONFIG.CONST_NUMERIC:
                  assert(
                    isinstance(var_info['consequences'][plugin],(int,float,complex))
                  )
                  if CONFIG.CONST_NUMERIC[plugin] == 'higher':
                    if var_info['consequences'][plugin] > plugin_criteria:
                      var_with_possible_annotations.append(mask_name)
                  elif CONFIG.CONST_NUMERIC[plugin]=='lower':
                    if var_info['consequences'][plugin] < plugin_criteria:
                      var_with_possible_annotations.append(mask_name)
                  else:
                    assert(False),f"error comparing numeric values for {gene} {var} {plugin} {study}"
                else:
                  assert(False),f"error values not numeric/categorical {gene} {var} {plugin} {study}"
        # assign variant to a mask based on order
        study_plugin_orders = CONFIG.plugin_orders[study]
        if len(var_with_possible_annotations) > 0:
          var_mask = None
          for i in var_with_possible_annotations:
            if i in study_plugin_orders:
              var_mask = i
              break
          assert(var_mask != None),f"did not assign unique mask {gene} {var} {study} {var_with_possible_annotations}"
          study_ofile = os.path.join(WDIR,study,f"{VCF_NAME}_annotations.txt")
          with open(study_ofile,'a') as ptr:
            ptr.write(
              f"{var}\t{gene}\t{var_mask}\n"
            )
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




  main()