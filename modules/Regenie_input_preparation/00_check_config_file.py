"""
Check the variables specified in the config file to make sure it has what it needs to run the scripts
"""
import yaml,os,shutil,sys,argparse
from collections import namedtuple
from itertools import compress
from pathlib import Path


#%% General
# check all the programs exists
def __check_path_exists():
  """check to make sure the required paths exists
  """
  path_existence = {
    'bcftools': os.path.isfile(CONFIG.bcftools), 
    "tabix": os.path.isfile(CONFIG.tabix), 
    "apptainer": os.path.isfile(CONFIG.apptainer),
    "plink" : os.path.isfile(CONFIG.plink),
    "plink2" : os.path.isfile(CONFIG.plink2),
    "vep_docker": os.path.isfile(CONFIG.vep_docker_image),
    "wdir": os.path.isdir(WDIR)
  }
  assert(
    all(list(path_existence.values()))
  ),f"These are not found: {';'.join(list(compress(list(path_existence.keys()),[not x for x in list(path_existence.values())])))}"

  return


def __check_build():
  """Check to make sure the build is one that we work with (only GRCh38)
  ...might have to add other things to verify build (like check positions and what not...)
  """
  assert(
    CONFIG.genome_build == 'GRCh38'
  ),f"specified {CONFIG.genome_build}, not supported"
  return


def  __check_numeric_constant():
  """make sure the numeric constant values are either higher or lower
  """
  for k,v in CONFIG.CONST_NUMERIC.items():
    assert(
      v in ['higher','lower']
    ),f"Not 'higher' or 'lower' for Numeric constant plugin {k}"
  return


#%%
def main():
  __check_path_exists()
  __check_build()
  __check_numeric_constant()

  print("Configuration file ok")

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

  if cargs.test == 't':
    from unittest import mock
    cargs = mock.Mock()
    cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"
    cargs.wdir="/scratch/richards/kevin.liang2/exwas_pipeline/results/pipeline_results"
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/sitesonly_VCF/wes_qc_chr16_sitesonly.vcf"
    print("TEST")


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(os.path.isfile(cargs.input_vcf)),'input vcf is missing'
  assert(cargs.wdir),'output directory missing'
  
  print("Checking configuration file")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"input VCF: {os.path.basename(cargs.input_vcf)}")
  print(f"output dir: {cargs.wdir}")
  print("="*20)


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  VCF_NAME = Path(cargs.input_vcf).stem
  WDIR = cargs.wdir

  main()