"""
  This is base on Ethan's script but rewritten in python with inputs specified from configuration file
"""
import os,yaml,pickle,pyreadr,shutil,sys
import pandas as pd
import numpy as np
from collections import namedtuple
import subprocess as sp
import argparse


def annotate_chr(vcf_infile,vcf_anno_out):
  # binds the cache directory
  apptainer_cmd = [
    CONFIG.apptainer,
    # set up apptainer runtime file system access
    # -C nothing in host system accessible except for binded ones
    # https://docs-dev.alliancecan.ca/wiki/Using_Apptainer#Official_Apptainer_documentation
    'run',"-C",'--bind',f"{CONFIG.vep_cache_dir}:/tmp/vep_cache,{WDIR}:/tmp/vep_wdir,{CONFIG.vep_plugin_dir}:/tmp/vep_plugins",
    # for input format
    'vep','-i',vcf_infile,
    '--assembly',CONFIG.genome_build,
    '--format','vcf', 
    '--cache',
    # corresonds to directories binded above
    '-o',f"/tmp/vep_wdir/{vcf_anno_out}",
    '--dir_cache','/tmp/vep_cache',
    '--dir_plugins','/tmp/vep_plugins'
  ]
  for x in [['--plugin',x] for x in CONFIG.vep_plugins]:
    apptainer_cmd += x
  apptainer_cmd += [
    '--offline',
    '--symbol','--coding_only','--no_stats'
    '--fork',"1"
  ]
  print("VEP annotation apptainer commands:")
  print(" ".join(apptainer_cmd))
  print("="*20)

  apptainer_run = sp.run(
    apptainer_cmd,
    stderr = sp.PIPE,
    check=True
  )

  print(apptainer_run.stderr.decode('utf-8'))
  print("="*20)

  return

def main():
  vcf_infile = os.path.join(WDIR,f'1_{VCF_NAME}_bcftool_variant_only_vcf.set_id.no_genotypes')
  vcf_anno_out=f'2_{VCF_NAME}_vcf_final_annotation.txt'
  
  annotate_chr(vcf_infile,vcf_anno_out)


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

  if cargs.test == 't':
    import mock
    cargs = mock.Mock()
    cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"
    cargs.wdir="/scratch/richards/kevin.liang2/exwas_pipeline/results/pipeline_results"
    cargs.input_vcf="/scratch/richards/kevin.liang2/exwas_pipeline/data/example_homo_sapiens_GRCh38.vcf"
    print("TEST")


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(os.path.isfile(cargs.input_vcf)),'input vcf is missing'
  assert(cargs.wdir),'output directory missing'

  print("Annotating VCF")
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
  
  main()

