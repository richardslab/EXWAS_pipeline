""" Create the set list file, which is shared by all studies. It is based on the input VCF.

"""

import shutil,os,sys,yaml,pyreadr,re,json,gzip,sqlite3
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse
from pathlib import Path

def main():
  print("creating set list file")
  print("*"*20)
  vep_summarie_file = os.path.join(
    WDIR,f"5_1_{VCF_NAME}_vep_summaries.sqlite3.db"
  )
  assert(
    os.path.isfile(vep_summarie_file)
  ),f"missing summary file"


  conn = sqlite3.connect(vep_summarie_file)
  cur = conn.cursor()
  genes = cur.execute(
    """
    SELECT distinct gene FROM vep_summaries
    """
  ).fetchall()

  setlist_file = os.path.join(
    WDIR,f"6_{VCF_NAME}.setlist"
  )
  print(f"set list file written to {setlist_file}")
  with open(setlist_file,'w') as ptr:  
    for gene in genes:
      gene = gene[0] # select returns a list of tuples
      var_var_location = cur.execute(
        """
        SELECT distinct SNP,location FROM vep_summaries WHERE gene = :gene
        """,{"gene":gene}
      ).fetchall()
      min_position = None
      gene_chr = None
      all_vars = []
      for var,var_location in var_var_location:
        all_vars.append(var)
        var_location_info = var_location.split(":")
        assert(
          len(var_location_info) == 2
        ),f"unknown SNP location notation {var_location}"
        
        # for indels, positions noted as start-end use the first position
        if (re.search("-",var_location_info[1])):
          var_pos = int(var_location_info[1].split("-")[0])
        else:  
          var_pos = int(var_location_info[1])

        var_chr = str(re.sub("chr","",var_location_info[0]))
        if min_position == None:
          min_position = var_pos
        else:
          if var_pos < min_position:
            min_position = var_pos
        if gene_chr == None:
          gene_chr = var_chr
        else:
          assert(gene_chr == var_chr),f"different chr for same gene {gene} {var}"
      all_vars_str = ",".join(all_vars)
      print(f"{gene}\t{gene_chr}\t{min_position}\t{all_vars_str}\n")
      ptr.write(
        f"{gene}\t{gene_chr}\t{min_position}\t{all_vars_str}\n"
      )
  conn.close()

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
    from unittest import mock
    cargs = mock.Mock()
    cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"
    cargs.wdir="/scratch/richards/kevin.liang2/exwas_pipeline/results/pipeline_results"
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/pipeline_results/example_homo_sapiens_GRCh38.vcf"
    print("TEST")


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  VCF_NAME = Path(cargs.input_vcf).stem
  WDIR = cargs.wdir


  print("Creating se list file")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"input VCF: {os.path.basename(cargs.input_vcf)}")
  print(f"output dir: {cargs.wdir}")
  print("="*20)





  main()