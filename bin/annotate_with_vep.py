#!/usr/bin/env python3
"""
  Annotate VCF files using VEP
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
    'run',"-C",'--bind',f"{CONFIG.vep_cache_dir}:/tmp/vep_cache,{IDIR}:/tmp/vep_wdir,{WDIR}:/tmp/vep_odir,{CONFIG.vep_plugin_dir}:/tmp/vep_plugins",
    VEP_IMG
  ]
  ## all the paths in here are relative to the paths within the apptainer image
  vep_cmd = [
    # for input format
    'vep','-i',f"/tmp/vep_wdir/{vcf_infile}",
    '--assembly',CONFIG.genome_build,
    '--format','vcf', 
    '--cache',
    "--force_overwrite",
    # corresonds to directories binded above
    '-o',f"/tmp/vep_odir/{vcf_anno_out}",
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
  vep_cmd += [
    '--offline',
    '--symbol',
    "--stats_text",
    "--compress_output",'bgzip'
  ]
  if CONFIG.vep_fork is not None:
    vep_cmd += [
      '--fork',str(CONFIG.vep_fork)
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
  print("Finished annotation")
  print("="*20)
  return

def main():
  vcf_infile = f'2_bcftool_sitesonly_{VCF_NAME}.vcf.gz'
  vcf_anno_out=f'3_annotation_results_{VCF_NAME}.txt.gz'
  
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
    '--input_dir',
    dest='input_dir',
    help="input directory of Sites only VCF",
    type=str
  )
  parser.add_argument(
    '--vep_img',
    dest='vep_img',
    help="VEP apptainer img",
    type=str
  )
  cargs =   parser.parse_args()

  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)

  VCF_NAME = Path(cargs.input_vcf).stem
  WDIR = os.getcwd()
  VEP_IMG = cargs.vep_img
  IDIR = cargs.input_dir
  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(os.path.isfile(cargs.input_vcf)),'input vcf is missing'

  print("Annotating VCF")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"input VCF: {os.path.basename(cargs.input_vcf)}")
  print("="*20)

  
  
  main()

