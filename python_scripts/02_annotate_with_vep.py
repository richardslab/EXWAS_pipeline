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
  apptainer_cmd = [
    CONFIG.apptainer,
    'run','--bind',f"{CONFIG.wdir}:{CONFIG.wdir}",
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
  vcf_infile = os.path.join(CONFIG.wdir,'1_bcftool_variant_only_vcf.set_id.no_genotypes')
  vcf_anno_out=os.path.join(CONFIG.wdir,'2_vcf_final_annotation.txt')
  
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
  cargs =   parser.parse_args()


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  print(f"Using {os.path.basename(cargs.cfile)}")


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  main()

