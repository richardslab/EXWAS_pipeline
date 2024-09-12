"""
  This is base on Ethan's script but rewritten in python with inputs specified from configuration file
"""
import os,yaml,pickle,pyreadr,shutil,sys
import pandas as pd
import numpy as np
from collections import namedtuple
import subprocess as sp
import argparse

def normalize_vcf(vcf_outfile):
  # commands to align, change the SNP IDs and remove sample info
  bcftool_align_cmd = [
    CONFIG.bcftools,
    'norm','-m','-any','--check-ref','w','-f',
    CONFIG.ref_alignment,
    CONFIG.input_vcf,
    '-Ou'
  ]
  bcftool_annotate_cmd = [
    CONFIG.bcftools,
    '--set-id','%CHROM:%POS:%REF:$FIRST_ALT',
    '-Ou'
  ]
  bcftool_var_only_file_cmd = [
    CONFIG.bcftools,
    'view','--drop-genotypes','-Oz'
  ]

  # run the commands and create the file
  bcftool_align = sp.Popen(
    bcftool_align_cmd,
    stdout=sp.PIPE,
    stderr=sp.PIPE,
    universal_newlines=True
  )
  bcftool_annotate = sp.Popen(
    bcftool_annotate_cmd,
    stdin = bcftool_align.stdout,
    stdout = sp.PIPE,
    universal_newlines=True
  )
  with open(vcf_outfile,'w') as ptr:
    bcftool_gen_var_only = sp.Popen(
      bcftool_var_only_file_cmd,
      stdin = bcftool_annotate.stdout,
      stdout = ptr,
      universal_newlines=True
    )
    # when each process finishes, sends signal to indicate it is done
    bcftool_align.stdout.close()
    bcftool_annotate.stdout.close()
    bcftool_gen_var_only.stdout.close()
    # wait for all processes to finish and check exit status
    bcftool_align.wait()
    bcftool_annotate.wait()
    bcftool_gen_var_only.wait()
    assert(
      bcftool_align.returncode == 0 and
      bcftool_annotate.returncode == 0 and
      bcftool_gen_var_only.returncode == 0
    ),'issue with generating vcf'

  # index the vcf file with tabix
  tabix_cmd = [
    CONFIG.tabix,
    '-p','vcf',vcf_outfile
  ]
  sp.run(
    tabix_cmd,check=True
  )

  return

def generate_plink_files(vcf_outfile,plink_output):
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
  vcf_outfile = os.path.join(WDIR,f'1_{VCF_NAME}_bcftool_variant_only_vcf.set_id.no_genotypes')
  normalize_vcf(
    vcf_outfile
  )
  # plink_output = os.path.join(WDIR,f'1_{VCF_NAME}_bcftool_variant_only_vcf.set_id.no_genotypes') 
  # generate_plink_files(vcf_outfile,plink_output)

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

