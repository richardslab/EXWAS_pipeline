#!/usr/bin/env python3
"""
  Find the paths and files for all the results

"""
import os,shutil,yaml,re,glob,itertools,gzip,argparse
import pandas as pd
from collections import namedtuple


def __find_files_per_phenotypes(each_study,unique_phenotypes):
  study_dir = os.path.join(WDIR,each_study)
  summary_out = os.path.join(study_dir,'Regenie_Summaries')
  os.makedirs(summary_out,exist_ok=True)
  regenie_results = glob.glob(
    os.path.join(REGENIE_S2_DIR,each_study,"8_regenie_S2_OUT_*.regenie.gz")
  )
  print(f"Processing {each_study}")
  print(f"results stored in {summary_out}")
  print("*"*20)
  files_per_pheno = {}
  for each_pheno in unique_phenotypes:
    phenotype_files = list(itertools.compress(
      regenie_results,
      [
      bool(re.search(f"_{each_pheno}\\.regenie\\.gz",x)) for x in regenie_results]
    ))
    if len(phenotype_files) == 0:
      assert(False),f"{each_study} missing files for {each_pheno}"
    files_per_pheno.update({each_pheno:phenotype_files})
  print(f"Found:")
  for k,v in files_per_pheno.items():
    print(f"phenotype: {k}. Number of files: {len(v)}")
  print("="*20)
  with gzip.open(os.path.join(summary_out,'Phenotypes_results_paths.yaml.gz'),'wt') as ptr:
    yaml.dump(files_per_pheno,ptr)
  return

def main():
  studies = list(CONFIG.mask_definitions.keys())
  unique_phenotypes = CONFIG.s2_params['--phenoCol'].split(",")
  print(f"Number of studies: {len(studies)}")
  print("\n".join([f"**{x}" for x in studies]))
  print(f"Number of phenotypes per study: {len(unique_phenotypes)}")
  print("\n".join([f"**{x}" for x in unique_phenotypes]))
  print("="*20)

  for each_study in studies:
    __find_files_per_phenotypes(each_study,unique_phenotypes)

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
    '--idir',
    dest='idir',
    type=str,
    help='Regenie S2 output directory'
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
  REGENIE_S2_DIR = cargs.idir

  main()