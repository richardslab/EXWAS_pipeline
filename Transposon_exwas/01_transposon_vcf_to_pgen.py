"""
Convert UKB VCF files to PGEN files
"""
import os,shutil,yaml,re,glob,csv,argparse
import subprocess as sp
from collections import namedtuple
from pathlib import Path


def main():
  plink2_cmd = [
    args.plink2,
    '--vcf',cargs.input_vcf,
    "--double-id",
    "--make-pgen",
    '--output-chr','chrM',
    "--memory","15000",
    '--vcf-half-call','m',
    "--set-all-var-ids",'@:#:\\$r:\\$a',
    "--new-id-max-allele-len",
    "--out",os.path.join(WDIR,f"{VCF_NAME}")
  ]
  res = sp.run(
    plink2_cmd,check=True
  )
  assert(res.returncode == 0)
  return


if __name__ == "__main__":
  from unittest import mock
  cargs = mock.Mock()

  cargs.input_vcf = "/scratch/richards/guillaume.butler-laporte/transposons/original_final_joint_call/ukb_final_autosomes.vcf.gz"

  cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config_transposons.yml"

  cargs.wdir = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/ukb_transposon_results"
  
  with open(cargs.cfile,'r') as ptr:
    params = yaml.safe_load(ptr)['proj_config']

  args = namedtuple("params",params.keys())(**params)

  WDIR = cargs.wdir
  VCF_NAME = Path(cargs.input_vcf).stem
  main()