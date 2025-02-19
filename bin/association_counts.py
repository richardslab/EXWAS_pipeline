#!/usr/bin/env python3
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

def __get_counts_per_study(study_regenie_result_paths,n_phenos,phenotype):
  """Obtain the counts for each phenotype in each study.
  
  Study-wide Bonferroni threshold is use assuming 20_000 genes.
    * Bonferroni threshold = 0.05/20_000/N masks/N models/N phenotypes

  Args:
      study_regenie_result_paths (str): path to results. Obtained from 00_find_data.py
      phenotype (str): phenotype name
  """
  files = study_regenie_result_paths[phenotype]
  # obtain all the exwas results into a dataframe
  # but keep only those with p-value < 0.05/20_000/n_phenos
  # # this is because it is the maximum p-value allowable for bonferroni assuming 20_000 genes tested across n phenotypes
  # # this would be bonferroni of 1 mask, 1 burden test, 1 allele frequency.
  # # smaller dataframe -> so faster processing time.
  max_p_value = 0.05/20_000/n_phenos
  relevant_exwas = pd.DataFrame()
  # keeps track of number of masks and models to do multiple testing correction
  # Notes:
  # # Alt would include allele frequency threshold as it is the Mask per allele frequencies (singletons, 0.001, etc)
  masks = set()
  models = set()
  for each_f in files:
    res = pd.read_csv(each_f,sep="\t",comment="#",quoting=csv.QUOTE_NONE)
    masks = masks.union(set(res.Alt.unique().tolist()))
    models = models.union(set(res.Model.unique().tolist()))
    assert(
      all(
        [x in res.columns for x in CONFIG.regenie_expected_columns]
      )
    ),f"Regenie results missing required columns.\nExpected: {','.join(CONFIG.regenie_expected_columns)}\nFound: {','.join(res.columns)}"
    res[['Pval',"Effect","LCI_Effect","UCI_Effect"]] = res[['Pval',"Effect","LCI_Effect","UCI_Effect"]].astype(float)
    relevant_exwas = pd.concat(
      [
        relevant_exwas,
        res.query("Pval <= @max_p_value")
      ],axis=0
    ).reset_index(drop=True)
  study_specific_bonferroni = max_p_value/len(masks)/len(models)/n_phenos
  study_counts = pd.DataFrame(
    columns = [
      "Phenotype","Mask","Model","Number of Exome-wide significant hits",
      "Top hit (Gene)",
      "Top hit (Gene) Effect size",
      "Top hit (Gene) P-value",
      "Study-wide Bonferroni threshold"
    ]
  )
  mask_model_pairs = [(mask,model) for mask in masks for model in models]
  assert(len(mask_model_pairs) == len(masks) * len(models))
  for mask,model in mask_model_pairs:
    mask_model_result = relevant_exwas.query(
      f"Alt == '{mask}' & Model == '{model}' & Pval <= {study_specific_bonferroni}"
    )
    n_exome_hit = len(mask_model_result)
    top_hits = mask_model_result.sort_values(['Pval'],ascending=True).head(n=1)
    if len(top_hits) == 0:
      top_hit_gene = None
      top_hit_pval = np.nan
      top_hit_beta = np.nan
    else:
      top_hit_gene = re.sub(f"\.{mask}","",top_hits.Name.iloc[0])
      top_hit_pval = top_hits.Pval.iloc[0]
      top_hit_beta = top_hits.Effect.iloc[0]
    study_counts = pd.concat(
      [
        study_counts,
        pd.DataFrame(
          {
            "Phenotype":phenotype,
            "Mask":mask,
            "Model":model,
            "Number of Exome-wide significant hits":len(mask_model_result),
            "Top hit (Gene)": top_hit_gene,
            "Top hit (Gene) Effect size": top_hit_beta,
            "Top hit (Gene) P-value": top_hit_pval,
            "Study-wide Bonferroni threshold": study_specific_bonferroni
          },index=[0]
        )
      ],axis=0
    ).reset_index(drop=True)
  return study_counts

def get_counts(each_study):
  """Obtain counts for each study.

  Phenotypes are processed in parallele.

  Args:
      each_study (str): study name specified in configuration file.
  """
  studies = list(CONFIG.mask_definitions.keys())
  unique_phenotypes = CONFIG.s2_params['--phenoCol'].split(",")
  assert(len(RES_PATH) == len(studies)),f"studies without results?"
  for each_study_file in RES_PATH:
    study_name = os.path.basename(os.path.normpath(re.sub(os.path.join("Regenie_Summaries","Phenotypes_results_paths.yaml.gz"),"",each_study_file)))
    print(f"Working on {each_study_file}")
    print(f"Working on {study_name}")
    assert(study_name in studies),f"unknown {study_name}" 
    with gzip.open(each_study_file,'rt') as ptr:
      study_regenie_result_paths = yaml.safe_load(ptr)
    summary_out = os.path.join(WDIR,each_study,'Regenie_Summaries')
    os.makedirs(summary_out,exist_ok=True)
    phenotypes = list(study_regenie_result_paths.keys())
    with Pool(PROCESSING_THREADS) as p:
      all_res = p.map(
        partial(__get_counts_per_study,study_regenie_result_paths,len(phenotypes)),
        phenotypes
      )
  result_counts  = pd.DataFrame()
  for res in all_res:
    columns = res.columns.tolist()
    res['Study'] = each_study
    res = res[["Study"]+columns]
    result_counts = pd.concat(
      [
        result_counts,
        res
      ],axis=0
    ).reset_index(drop=True)
  # format the display
  numeric_cols = [
    "Top hit (Gene) Effect size",
    "Top hit (Gene) P-value",
    "Study-wide Bonferroni threshold"
  ]
  for numeric_col in numeric_cols:
    round_digits = (np.abs(result_counts[numeric_col]) >= 0.0001).astype(bool) & (~result_counts[numeric_col].isna()).tolist() & (np.abs(result_counts[numeric_col]) < 1).astype(bool)
    sci_digits = (result_counts[numeric_col] < 0.0001).astype(bool) & ~result_counts[numeric_col].isna() & (np.abs(result_counts[numeric_col]) > 0 ).astype(bool)
    result_counts.loc[round_digits,numeric_col] = np.round(result_counts.loc[round_digits,numeric_col],3)
    result_counts.loc[sci_digits,numeric_col] = result_counts.loc[sci_digits,numeric_col].apply(lambda val: "{:.3e}".format(val))
  with gzip.open(os.path.join(summary_out,'ExWAS_counts.tsv.gz'),'wt') as ptr: 
    result_counts.to_csv(ptr,quoting=csv.QUOTE_NONE,header=True,index=False,sep="\t",na_rep="NA")
  return 

def main():
  studies = list(CONFIG.mask_definitions.keys())
  for each_study in studies:
    get_counts(each_study)
  



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
    __file__ = "/home/richards/kevin.liang2/scratch/exwas_pipeline/src/modules/Regenie_results_summaries/01_compute_lambda.py"
    print("TEST")


  assert(os.path.isfile(cargs.cfile)),'config file is missing'

  print("Obtain summary counts")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print("Study wide Bonferroni corrected threshold assuming 20K genes are used")
  print("**i.e., Bonferroni threshold = 0.05/20K/number of masks/number of models for each phenotype")
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