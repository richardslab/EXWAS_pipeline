#!/usr/bin/env python3
""" 
Create the VEP apptainer image
"""

import shutil,os,sys,yaml,re,json,gzip,sqlite3
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse
from pathlib import Path
import subprocess as sp

def main():
  build_cmd = [
    CONFIG.apptainer,'build',
    "--fakeroot",
    "--build-arg",f"vep_version={str(CONFIG.vep_vrs)}",
    "vep_apptainer.sif",apptainer_def
  ]
  _ = sp.run(
    build_cmd,check=True
  )
  return 

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--config_file','-c',
    dest='cfile',
    help='configuration yaml file'
  )
  parser.add_argument(
    '--vep_apptainer_def_file',
    dest='vep_def_file',
    help="VEP apptainer definition file",
    type=str
  )
  cargs =   parser.parse_args()


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  apptainer_def = cargs.vep_def_file
  WDIR = os.getcwd()


  print("Creating apptainer img")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  main()