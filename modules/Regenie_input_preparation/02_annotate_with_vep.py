"""
  This is base on Ethan's script but rewritten in python with inputs specified from configuration file
"""
import os,yaml,pickle,pyreadr,shutil,sys
import pandas as pd
import numpy as np
from collections import namedtuple
import subprocess as sp
import argparse


def annotate_chr(chr,vcf_infile,vcf_anno_out):
  # binds the cache directory
  apptainer_cmd = [
    CONFIG.apptainer,
    'run',"-C",'--bind',f"{CONFIG.vep_cache_dir}:/tmp/vep_cache,{WDIR}:/tmp/vep_wdir",
    'vep','-i',vcf_infile,
    '--assembly',CONFIG.genome_build,
    '--format','vcf',
    '--cache',
    '-o',vcf_anno_out,
    '--dir_cache',CONFIG.vep_cache
  ]
  for x in [['--plugin',x] for x in CONFIG.vep_plugins]:
    apptainer_cmd += x
  apptainer_cmd += [
    '--force_overwrite','--offline',
    '--symbol','--coding_only','--no_stats'
    '--fork',1,
    '--quiet'
  ]
  sp.run(apptainer_cmd,check=True)

  return

def main():
  vcf_infile = os.path.join(WDIR,f'1_{VCF_NAME}_bcftool_variant_only_vcf.set_id.no_genotypes')
  vcf_anno_out=os.path.join(WDIR,f'2_{VCF_NAME}_vcf_final_annotation.txt')
  
  annotate_chr(vcf_infile,vcf_anno_out)


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
  
  main()

