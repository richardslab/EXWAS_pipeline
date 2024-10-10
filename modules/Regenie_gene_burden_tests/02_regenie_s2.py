"""Run Step 2 of regenie using genotyped variants
"""
import os,sys,shutil,yaml,re,argparse
from collections import namedtuple
import subprocess as sp


def _run_regenie_s2_each_study(study):
  study_outdir = os.path.join(WDIR,study)
  assert(os.path.isdir(study_outdir)),f"missing a step for {study_outdir}??"

  s2_cmd = [
    CONFIG.regenie,
    "--step","2",
    "--anno-file",os.path.join(study_outdir,f"{VCF_NAME}_annotations.txt"),
    "--mask-def",os.path.join(study_outdir,f"{VCF_NAME}_masks.txt"),
    "--set-list",os.path.join(WDIR,f"6_{VCF_NAME}.setlist")
  ]
  for k,v in CONFIG.s2_params.items():
    if k == "--lowmem-prefix":
      v = os.path.join(WDIR,v)
    if v == "":
      s2_cmd += [k]
    elif v != "":
      s2_cmd += [f"{k} {v}"]
    else:
      assert False, 'invalid Regenie Step 2 value format issue in config file {k} {v}'    
  s2_cmd += [
    "--pred",os.path.join(WDIR,f"7_{VCF_NAME}_regenie_S1_OUT_pred.list"),
    "--out",os.path.join(study_outdir,f"8_{VCF_NAME}_regenie_S2_OUT")
  ]
  print("Regenie S2 command:")
  print(" ".join(s2_cmd))
  print("*"*20)
  s2_out = sp.run(
    " ".join(s2_cmd),check=True,shell=True,stderr=sp.PIPE
  )
  print(s2_out.stderr.decode("utf-8"))
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
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/data/wes_qc_chr3_chr_full_final.vcf.subset.sorted.vcf.gz"
    print("TEST")
  



  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(os.path.isfile(cargs.input_vcf)),'input vcf is missing'
  assert(cargs.wdir),'output directory missing'
  
  print("Run Regenie Step 2")
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