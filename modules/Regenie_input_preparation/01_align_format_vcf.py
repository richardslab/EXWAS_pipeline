"""
  This is base on Ethan's script but rewritten in python with inputs specified from configuration file
"""
import os,yaml,pickle,pyreadr,shutil,sys,gzip
import pandas as pd
import numpy as np
from collections import namedtuple
import subprocess as sp
import argparse


def normalize_vcf(vcf_outfile):
  """Normalize input VCF against reference assembly

  Args:
      vcf_outfile (str): input vcf file
  """
  # commands to align, change the SNP IDs and remove sample info
  bcftool_align_cmd = [
    CONFIG.bcftools,
    'norm','-m','-any','--check-ref','w','-f',
    CONFIG.reference_fasta,
    cargs.input_vcf,
    '-O','u'
  ]
  bcftool_annotate_cmd = [
    CONFIG.bcftools,
    'annotate',
    '--set-id','%CHROM:%POS:%REF:%FIRST_ALT',
    '-O','u'
  ]
  bcftool_var_only_file_cmd = [
    CONFIG.bcftools,
    'view','--drop-genotypes','-O','z','-o',vcf_outfile
  ]
  # print the equivalent bash command to log
  print("Alignment Command: ")
  print(
    f"{' '.join(bcftool_align_cmd)} | {' '.join(bcftool_annotate_cmd)} | {' '.join(bcftool_var_only_file_cmd)}"
  )
  print("="*20)
  # run the commands and create the file
  bcftool_align = sp.Popen(
      bcftool_align_cmd,
      stdout = sp.PIPE,
      stderr = sp.PIPE,
      universal_newlines=True
    )
  bcftool_annotate = sp.Popen(
    bcftool_annotate_cmd,
    stdin = bcftool_align.stdout,
    stdout = sp.PIPE,
    universal_newlines=True
  )
  bcftool_gen_var_only = sp.run(
    bcftool_var_only_file_cmd,
    stdin = bcftool_annotate.stdout,
    check=True
  )
  # when each process finishes, sends signal indicate it is done
  bcftool_align.stdout.close()
  bcftool_annotate.stdout.close()
  # wait for all processes to finish and check exit status
  bcftool_align_res = bcftool_align.communicate()
  bcftool_annotate_res = bcftool_annotate.communicate()
  assert(
    bcftool_align.returncode == 0 and
    bcftool_annotate.returncode == 0 and
    bcftool_gen_var_only.returncode == 0
  ),'issue with generating vcf'
  print(bcftool_align_res.stdout.decode('utf-8'))
  print(bcftool_align_res.stderr.decode('utf-8'))
  print(bcftool_annotate_res.stdout.decode('utf-8'))
  print(bcftool_annotate_res.stderr.decode('utf-8'))
  print(bcftool_gen_var_only.stdout.decode('utf-8'))
  print(bcftool_gen_var_only.stderr.decode('utf-8'))
  

  # index the vcf file with tabix
  tabix_cmd = [
    CONFIG.tabix,
    '-p','vcf',vcf_outfile
  ]
  print("Tabix Command: ")
  print(
    " ".join(tabix_cmd)
  )
  print("="*20)
  tabix_run = sp.run(
    tabix_cmd,check=True
  )
  print(tabix_run.stdout.decode("utf-8"))
  print(tabix_run.stderr.decode("utf-8"))
  
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
  vcf_outfile = os.path.join(WDIR,f'1_{VCF_NAME}_bcftool_variant_only.set_ids.no_genotypes.vcf.gz')
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
    cargs.wdir="/scratch/richards/kevin.liang2/exwas_pipeline/results/pipeline_results"
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/data/wes_qc_chr3_chr_full_final.vcf.subset.sorted.vcf.gz"
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
  VCF_NAME = os.path.basename(cargs.input_vcf)
  WDIR = cargs.wdir

  main()

