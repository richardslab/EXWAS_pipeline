"""
  Generate some summaries from the ExWAS

  1. Computes the multiple testing burden per study
    * this is based on 0.5/20K genes/N masks/N models/N phenotypes

  - The number of significant associations
  - The strongest association for each trait under each condition
"""

import os,shutil,yaml,re,glob,itertools,gzip,argparse,csv,sys
import pandas as pd
import numpy as np
from collections import namedtuple
from multiprocessing import Pool,cpu_count
from functools import partial

def get_counts():
  study_dir = os.path.join(WDIR,each_study)
  regenie_s2_dir = os.path.join(study_dir,'Regenie_S2')
  summary_out = os.path.join(study_dir,'Regenie_Summaries')
  summary_file = os.path.join(summary_out,'Phenotypes_results_paths.yaml.gz')
  assert(os.path.isfile(summary_file)),f"Missing phenotype paths. Run 00_find_data.py or error in earlier steps"
  with gzip.open(summary_file,'rt') as ptr:
    study_regenie_result_paths = yaml.safe_load(ptr)

  return
def main():
  studies = list(CONFIG.mask_definitions.keys())
  unique_phenotypes = CONFIG.s2_params['--phenoCol'].split(",")


  return

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--config_file','-c',
    dest='cfile',
    type=str,
    help='configuration yaml file'
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
    import mock
    cargs = mock.Mock()
    cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"
    cargs.wdir="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/PAST_pipeline_results/intial_runs"
    __file__ = "/home/richards/kevin.liang2/scratch/exwas_pipeline/src/modules/Regenie_results_summaries/01_compute_lambda.py"
    print("TEST")


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(cargs.wdir),'output directory missing'

  print("Obtain summary counts")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"output dir: {cargs.wdir}")
  print("="*20)


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  WDIR = cargs.wdir
  if not isinstance(CONFIG.processing_threads,(float,int)):
    PROCESSING_THREADS = np.min([5,cpu_count()-1])
  else:
    PROCESSING_THREADS = np.min([CONFIG.processing_threads,cpu_count()-1])

  sys.path.append(os.path.dirname(__file__))
  from python_helpers import exwas_helpers

  main()