#!/usr/bin/env python3
"""Run Step 2 of regenie using genotyped variants
"""
import os,sys,shutil,yaml,re,argparse
from collections import namedtuple
import subprocess as sp
import warnings
from pathlib import Path

def _run_regenie_s2_each_study(study):
  regenie_s1_dir = os.path.dirname(REGENIE_S1_RES)
  study_outdir = os.path.join(WDIR,study)
  input_file_dir = os.path.join(IDIR,study)
  assert(os.path.isdir(input_file_dir)),f"missing a step for {input_file_dir}??"
  regenie_s2_dir = os.path.join(study_outdir)
  os.makedirs(regenie_s2_dir,exist_ok=True)

  print("Looking for relevant supplementary files")
  expected_annotation_file,expected_mask_file,expected_setlist_file,wildcard_characters = regenie_helpers.__find_regenie_supplementary_files(
    input_vcf = cargs.input_vcf,
    nxtflow_g = cargs.nxtflow_g,
    nxtflow_annotation = cargs.nxtflow_annotation,
    study_dir = input_file_dir
  )
  print("*"*20)
  print("found the following files")
  print(f"Regenie Input file: {cargs.input_vcf}.{cargs.nxtflow_gtype}")
  print(f"Annotation file: {expected_annotation_file}")
  print(f"Set list file: {expected_setlist_file}")
  print(f"Mask file: {expected_mask_file}")
  print("="*20)

  s2_cmd = [
    CONFIG.regenie,
    "--step","2",
    "--anno-file",expected_annotation_file,
    "--mask-def",expected_mask_file,
    "--set-list",expected_setlist_file
  ]
  if cargs.nxtflow_gtype == "pgen":
    s2_cmd += [
      "--pgen",cargs.input_vcf
    ]
  elif cargs.nxtflow_gtype == "bgen":
    s2_cmd += [
      "--bgen",f"{cargs.input_vcf}.bgen"
    ]
    if "--sample" not in CONFIG.s2_params:
      # warnings goes to stderr, so print also to put it in log
      # can also redirect outputs in pipeline
      warnings.warn("Using bgen files without sample files")
      print("WARNING: Using bgen files without sample files")
  elif cargs.nxtflow_gtype == "bed":
    c2_cmd += [
      "--bed",cargs.input_vcf
    ]
  else:
    raise ValueError("unrecognized genetic file type")
  for k,v in CONFIG.s2_params.items():
    if k == "--lowmem-prefix":
      v = os.path.join(WDIR,v)
    if v == "":
      s2_cmd += [k]
    elif v != "":
      s2_cmd += [k,v]
    else:
      assert False, 'invalid Regenie Step 2 value format issue in config file {k} {v}'    
  # add the --htp flag to get case counts and standardized output format
  if "--htp" not in CONFIG.s2_params:
    s2_cmd += [
      "--htp","htp_results"
    ]
  s2_cmd += [
    "--pred",REGENIE_S1_RES,
    "--out",os.path.join(regenie_s2_dir,f"8_regenie_S2_OUT_{VCF_NAME}")
  ]
  print("Regenie S2 command:")
  print(" ".join(s2_cmd))
  print("*"*20)
  s2_out = sp.run(
    s2_cmd,check=True,stderr=sys.stderr
  )
  print("="*20)

  return

def main():
  studies = list(CONFIG.mask_definitions.keys())
  assert(len(list(set(studies))) == len(studies)),"duplicated studies specified"

  for each_study in studies:
    _run_regenie_s2_each_study(each_study)

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
    '--regenie_s1',
    dest='regenie_s1',
    help="Regenie S1 pred list file full path",
    type=str
  )
  parser.add_argument(
    '--idir',
    dest='idir',
    help="Regenie input dir",
    type=str
  )
  parser.add_argument(
    "--nxtflow_genetic",
    dest="nxtflow_g",
    help="Nextflow 'step2_exwas_genetic' input parameter"
  )
  parser.add_argument(
    "--nxtflow_genetic_type",
    dest="nxtflow_gtype",
    help="Nextflow 'step2_exwas_genetic' input parameter"
  )
  parser.add_argument(
    "--nxtflow_annotation",
    dest="nxtflow_annotation",
    help="Nextflow 'annotation_vcf' input parameter"
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
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/data/exwas_data/wes_qc_chr10.bgen"
    cargs.nxtflow_gtype="bgen"
    cargs.nxtflow_g="/home/richards/kevin.liang2/scratch/exwas_pipeline/data/exwas_data/wes_qc_chr*"
    cargs.nxtflow_annotation = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/sitesonly_VCF/wes_qc_chr*_sitesonly.vcf"
    __file__ = "/home/richards/kevin.liang2/scratch/exwas_pipeline/src/modules/Regenie_gene_burden_tests/02_regenie_s2.py"
    print("TEST")
  



  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(os.path.isfile(cargs.input_vcf)),'input vcf is missing'
  
  print("Run Regenie Step 2")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"input VCF: {os.path.basename(cargs.input_vcf)}")
  print("="*20)


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  cargs.input_vcf = str(Path(cargs.input_vcf).with_suffix(""))
  VCF_NAME = Path(cargs.input_vcf).stem
  WDIR = os.getcwd()
  REGENIE_S1_RES = cargs.regenie_s1
  IDIR = cargs.idir

  sys.path.append(os.path.dirname(__file__))
  from python_helpers import regenie_helpers

  main()