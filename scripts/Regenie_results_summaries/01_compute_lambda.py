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

def __compute_lambda(each_study):
  study_regenie_result_paths,summary_out = exwas_helpers.__get_study_exwas_paths(WDIR,each_study)
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
  for each_study in studies:
    print(f"Working on {each_study}")
    __compute_lambda(each_study)

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
    cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/zhao_etal_config/proj_config.yml"
    cargs.wdir="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/zhao_etal_BSN_BMI"
    __file__ = "/home/richards/kevin.liang2/scratch/exwas_pipeline/src/modules/Regenie_results_summaries/01_compute_lambda.py"
    print("TEST")


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(cargs.wdir),'output directory missing'

  print("Computing genomic inflation factors (lambda)")
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