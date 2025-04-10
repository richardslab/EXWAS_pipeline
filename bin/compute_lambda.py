#!/usr/bin/env python3
"""
  Compute the lambda for each results
"""
import os,shutil,yaml,re,glob,itertools,gzip,argparse,csv,sys
import pandas as pd
import numpy as np
from collections import namedtuple
from multiprocessing import Pool,cpu_count
from functools import partial

def __process_single_phenotype_lambda(study_regenie_result_paths,phenotype):
  """Compute genomic inflation factor for each phenotype

  Will compute the genomic inflation factor across all chromosomes per phenotype per masks per model.

  Expect output columns from Regenie --htp flag.

  Genomic inflation factor computed as 
    median(observed chisq stat)/(median chisq stat with df = 1)
      * ratio of median of the test statistics
    
  Args:
      study_regenie_result_paths (str): path to the regenie results. Obtained from 00_find_data.py
      phenotype (str): Phenotype name

  """
  files = study_regenie_result_paths[phenotype]
  tmp_lambda_p_values = {}
  for each_file in files:
    file_cont = pd.read_csv(each_file,sep="\t",comment="#",quoting=csv.QUOTE_NONE)
    assert(
      all(
        [x in file_cont.columns for x in CONFIG.regenie_expected_columns]
      )
    ),f"Regenie results missing required columns.\nExpected: {','.join(CONFIG.regenie_expected_columns)}\nFound: {','.join(file_cont.columns)}"
    rel_content = file_cont[["Name","Trait","Alt","Model","Pval"]]
    unique_masks = rel_content.Alt.unique().tolist()
    for each_mask in unique_masks:
      mask_info = rel_content.query(f"Alt == '{each_mask}'")
      unique_models = mask_info.Model.unique().tolist()
      for each_model in unique_models:
        mask_model_pvalue = mask_info.query(
          f"Model == '{each_model}'"
        ).Pval.to_list()
        k = (each_mask, each_model)
        if k in tmp_lambda_p_values:
          tmp_lambda_p_values[k] += mask_model_pvalue
        else:
          tmp_lambda_p_values[k] = mask_model_pvalue
  regenie_lambdas = {}
  for each_k,pval in tmp_lambda_p_values.items():
    lambdaa = exwas_helpers.compute_lambda(pval)
    mask,model = each_k
    k = (phenotype,mask,model)
    regenie_lambdas[k] = lambdaa
  return regenie_lambdas

def __compute_lambda(each_study,each_study_file):
  with gzip.open(each_study_file,'rt') as ptr:
    study_regenie_result_paths = yaml.safe_load(ptr)
  summary_out = os.path.join(WDIR,each_study,'Regenie_Summaries')
  os.makedirs(summary_out,exist_ok=True)
  regenie_lambdas = {}
  phenotypes = list(study_regenie_result_paths.keys())
  with Pool(min(cpu_count(),PROCESSING_THREADS)) as p:
    all_regenie_results = p.map(
      partial(
        __process_single_phenotype_lambda,study_regenie_result_paths
      ),phenotypes
    )
  lambda_df = pd.DataFrame()
  for res in all_regenie_results:
    for k,v in res.items():
      phenotype,mask,model = k
      lambda_df = pd.concat(
        [
          lambda_df,
          pd.DataFrame(
            {
              "Phenotype":phenotype,
              "Mask":mask,
              "Model":model,
              "Lambda":v
            },index=[0]
          )
        ],axis=0
      ).reset_index(drop=True)
  with gzip.open(os.path.join(summary_out,"Genomic_inflation_factors.tsv.gz"),'wt') as ptr:
    lambda_df.to_csv(ptr,sep="\t",index=False,header=True,quoting=csv.QUOTE_NONE)
  return


def main():
  studies = list(CONFIG.mask_definitions.keys())
  unique_phenotypes = CONFIG.s2_params['--phenoCol'].split(",")
  assert(len(RES_PATH) == len(studies)),f"studies without results?"
  for each_study_file in RES_PATH:
    study_name = os.path.basename(os.path.normpath(re.sub(os.path.join("Regenie_Summaries","Phenotypes_results_paths.yaml.gz"),"",each_study_file)))
    print(f"Working on {each_study_file}")
    print(f"Working on {study_name}")
    assert(study_name in studies),f"unknown {study_name}" 
    __compute_lambda(study_name,each_study_file)

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
    '--res_path',
    dest='res_path',
    nargs="+",type=str,
    help='Path to Regenie results. Obtained from find_data.py'
  )
  cargs =   parser.parse_args()

  assert(os.path.isfile(cargs.cfile)),'config file is missing'

  print("Computing genomic inflation factors (lambda)")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print("="*20)


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  WDIR = os.getcwd()
  RES_PATH = cargs.res_path

  if not isinstance(CONFIG.processing_threads,(float,int)):
    PROCESSING_THREADS = np.min([5,cpu_count()-1])
  else:
    PROCESSING_THREADS = np.min([CONFIG.processing_threads,cpu_count()-1])

  sys.path.append(os.path.dirname(__file__))
  from python_helpers import exwas_helpers

  main()