#!/usr/bin/env python3
"""Run Step 1 of regenie using genotyped variants
"""
import os,sys,shutil,yaml,re,argparse
from collections import namedtuple
import subprocess as sp
from pathlib import Path

def main():
  regenie_s1_out = WDIR
  os.makedirs(regenie_s1_out,exist_ok=True)
  s1_cmd = [
    CONFIG.regenie,"--step","1"
  ]
  for k,v in CONFIG.s1_params.items():
    if k == "--lowmem-prefix":
      v = os.path.join(regenie_s1_out,v)
    if v == "":
      s1_cmd += [k]
    elif v != "":
      s1_cmd += [k,v]
    else:
      assert False, 'invalid Regenie step 1 value format issue in config file {k} {v}'
  s1_cmd += [
    "--out",os.path.join(regenie_s1_out,f"7_Regenie_S1")
  ]
  
  print("Regenie S1 command:")
  print(" ".join(s1_cmd))
  print("*"*20)

  s1_out = sp.run(
    s1_cmd,check=True,stderr= sys.stderr
  )
  print("="*20)
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
    '--test',
    default='f',
    type=str
  )
  cargs =   parser.parse_args()

  if cargs.test == 't':
    from unittest import mock
    cargs = mock.Mock()
    cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"
    print("TEST")
  



  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  
  print("Run Regenie Step 1")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print("="*20)


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  WDIR = os.getcwd()
  
  main()