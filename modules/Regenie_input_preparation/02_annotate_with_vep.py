"""
  This is base on Ethan's script but rewritten in python with inputs specified from configuration file
"""
import os,yaml,pickle,pyreadr,shutil,sys,re
import pandas as pd
import numpy as np
from collections import namedtuple
import subprocess as sp
import argparse
from pathlib import Path

def annotate_vcf(vcf_infile,vcf_anno_out):
  # binds the cache directory
  apptainer_cmd = [
    CONFIG.apptainer,
    # set up apptainer runtime file system access
    # -C nothing in host system accessible except for binded ones
    # https://docs-dev.alliancecan.ca/wiki/Using_Apptainer#Official_Apptainer_documentation
    'run',"-C",'--bind',f"{CONFIG.vep_cache_dir}:/tmp/vep_cache,{WDIR}:/tmp/vep_wdir,{CONFIG.vep_plugin_dir}:/tmp/vep_plugins",
    CONFIG.vep_docker_image
  ]

  ## all the paths in here are relative to the paths within the apptainer image
  vep_cmd = [
    # for input format
    'vep','-i',f"/tmp/vep_wdir/{vcf_infile}",
    '--assembly',CONFIG.genome_build,
    '--format','vcf', 
    '--cache',
    # corresonds to directories binded above
    '-o',f"/tmp/vep_wdir/{vcf_anno_out}",
    '--dir_cache','/tmp/vep_cache',
    '--dir_plugins','/tmp/vep_plugins'
  ]
  # add each plugin option and replace 
  #   $vep_cache_dir with /tmp/vep_cache
  #   $vep_plugins_dir with /tmp/vep_plugins
  for x in CONFIG.vep_plugins:
    replaced_x = re.sub("\$vep_cache_dir",'/tmp/vep_cache',x)
    replaced_x = re.sub("\$vep_plugins_dir",'/tmp/vep_plugins',replaced_x)
    vep_cmd += ['--plugin',replaced_x]
  # Do not use --fork option due to error with loftee
  # https://github.com/konradjk/loftee/issues/45
  vep_cmd += [
    '--offline',
    '--symbol','--coding_only',
    '--no_stats'
  ]

  full_cmd = apptainer_cmd + [f"\"{' '.join(vep_cmd)}\""]
  print("VEP annotation apptainer commands:")
  print(" ".join(full_cmd))
  print("*"*20)

  apptainer_run = sp.run(
    apptainer_cmd + [' '.join(vep_cmd)],
    stderr = sys.stderr,
    check=True
  )
  assert(apptainer_run.returncode == 0),f"Error with VEP"
  print("="*20)
  
  return

def main():
  vcf_infile = f'2_bcftool_sitesonly_{VCF_NAME}.vcf.gz'
  vcf_anno_out=f'3_annotation_results_{VCF_NAME}.txt'
  
  annotate_vcf(vcf_infile,vcf_anno_out)


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
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/sitesonly_VCF/wes_qc_chr16_sitesonly.vcf"
    print("TEST")

  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)

  VCF_NAME = Path(cargs.input_vcf).stem
  WDIR = cargs.wdir
  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(os.path.isfile(cargs.input_vcf)),'input vcf is missing'
  assert(cargs.wdir),'output directory missing'

  print("Annotating VCF")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"input VCF: {os.path.basename(cargs.input_vcf)}")
  print(f"output dir: {cargs.wdir}")
  print("="*20)

  
  
  main()

