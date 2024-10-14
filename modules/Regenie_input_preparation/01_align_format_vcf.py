"""
  This is base on Ethan's script but rewritten in python with inputs specified from configuration file
"""
import os,yaml,pickle,pyreadr,shutil,sys,gzip
import pandas as pd
import numpy as np
from collections import namedtuple
import subprocess as sp
import argparse
from pathlib import Path


def normalize_vcf(vcf_outfile):
  """Normalize input VCF against reference assembly

  Args:
      vcf_outfile (str): input vcf file
  """
  # commands to align, change the SNP IDs and remove sample info
  bcftool_var_only_file_cmd = [
    CONFIG.bcftools,
    'view','--drop-genotypes',
    '-Ou',cargs.input_vcf
  ]
  bcftool_align_cmd = [
    CONFIG.bcftools,
    'norm','-m','-any','--check-ref','w','-f',
    CONFIG.reference_fasta,
     '-Ou'
  ]
  bcftool_annotate_cmd = [
    CONFIG.bcftools,
    'annotate',
    '--set-id','%CHROM:%POS:%REF:%FIRST_ALT',
    '-Oz','-o',vcf_outfile
  ]
  # print the equivalent bash command to log
  print("Alignment Command: ")
  print(
    f"{' '.join(bcftool_var_only_file_cmd)} | {' '.join(bcftool_align_cmd)} | {' '.join(bcftool_annotate_cmd)}"
  )
  print("*"*20)

  # run the commands and create the file
  bcftool_gen_var_only = sp.Popen(
    bcftool_var_only_file_cmd,
    stderr = sp.PIPE,
    stdout = sp.PIPE
  )
  
  bcftool_align = sp.Popen(
      bcftool_align_cmd,
      stdout = sp.PIPE,
      stderr = sp.PIPE,
      stdin=bcftool_gen_var_only.stdout
    )

  bcftool_annotate = sp.Popen(
    bcftool_annotate_cmd,
    stdin = bcftool_align.stdout,
    stdout = sp.PIPE,
    stderr = sp.PIPE
  )


  bcftool_gen_var_only.stdout.close()
  bcftool_align.stdout.close()
  
  print(bcftool_gen_var_only.stderr.readlines())
  print(bcftool_align.stderr.readlines())
  print(bcftool_annotate.stderr.readlines())

  bcftool_gen_var_only.poll()
  bcftool_align.poll()
  bcftool_annotate.poll()

  bcftool_gen_var_only.wait()
  bcftool_align.wait()
  bcftool_annotate.wait()

  assert(
    bcftool_gen_var_only.returncode == 0 and
    bcftool_align.returncode == 0 and
    bcftool_annotate.returncode == 0
  ),'issue with generating vcf'
  
  print("="*20)

  # index the vcf file with tabix
  tabix_cmd = [
    CONFIG.tabix,
    '-p','vcf',vcf_outfile
  ]
  print("Tabix Command: ")
  print(
    " ".join(tabix_cmd)
  )
  print("*"*20)
  tabix_run = sp.run(
    tabix_cmd,check=True,stderr=sp.PIPE
  )
  print(tabix_run.stderr.decode("utf-8"))
  print("="*20)
  
  return

def generate_plink_files(vcf_outfile,plink_output):
  """Convert the vcf files into plink files

  Args:
      vcf_outfile (str): VCF file path 
      plink_output (str): plink output file paths
  """
  plink_cmd = [
    CONFIG.plink,
    '--vcf','--double-id',
    vcf_outfile,'--make-bed','--out',plink_output
  ]
  plink_prune_cmd = [
    CONFIG.plink,
    '--bfile',plink_output,
    '--hwe','1E-15 midp',
    '--maf','0.01','--geno','0.1',
    '--indep-pairwise','50 5 0.5',
    '--out','pruned_variants.txt'
  ]
  sp.run(
    plink_cmd,check=True
  )
  sp.run(plink_prune_cmd,check=True)

  return

def main():
  vcf_outfile = os.path.join(WDIR,f'2_bcftool_sitesonly_{VCF_NAME}.vcf.gz')
  normalize_vcf(vcf_outfile)

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
  
  if cargs.test =='t':
    import mock
    cargs = mock.Mock()
    cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"
    cargs.wdir="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/pipeline_results"
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/sitesonly_VCF/wes_qc_chr1_sitesonly.vcf"
    print("TEST")


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(os.path.isfile(cargs.input_vcf)),'input vcf is missing'
  assert(cargs.wdir),'output directory missing'

  print("aligning VCF files")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"input VCF: {os.path.basename(cargs.input_vcf)}")
  print(f"output dir: {cargs.wdir}")
  print("="*20)


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  VCF_NAME = Path(cargs.input_vcf).stem
  WDIR = cargs.wdir

  main()

