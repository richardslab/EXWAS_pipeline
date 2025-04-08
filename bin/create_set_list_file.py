#!/usr/bin/env python3
""" Create the set list file, which is shared by all studies. It is based on the input VCF.

"""

import shutil,os,sys,yaml,pyreadr,re,json,gzip,sqlite3
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse
from pathlib import Path
from itertools import compress

def main():
  ANNOTATION_FILE = os.path.join(IDIR,f"5_1_vep_summaries_{VCF_NAME}.sqlite3.db")
  assert(
    os.path.isfile(ANNOTATION_FILE)
  ),f"missing summary file {ANNOTATION_FILE}"


  try:
    conn = sqlite3.connect(ANNOTATION_FILE)
    cur = conn.cursor()
    genes = cur.execute(
      """
      SELECT distinct gene FROM vep_summaries
      """
    ).fetchall()
    conn.close()
  except Exception as e:
    print(f"SQLITE3 error: {e}")
    conn.close()
    raise

  setlist_file = os.path.join(
    WDIR,f"6_{VCF_NAME}.setlist"
  )
  print(f"set list file written to {setlist_file}")
  with open(setlist_file,'w') as ptr:  
    for gene in genes:
      gene = gene[0] # select returns a list of tuples
      try:
        conn = sqlite3.connect(ANNOTATION_FILE)
        cur = conn.cursor()
        var_var_location = cur.execute(
          """
          SELECT distinct SNP,location FROM vep_summaries WHERE gene = :gene
          """,{"gene":gene}
        ).fetchall()
        conn.close()
      except Excpetion as e:
        print(f"SQLITE3 error: {e}")
        conn.close()
        raise
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
      ptr.write(
        f"{gene}\t{gene_chr}\t{min_position}\t{all_vars_str}\n"
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
    '--input_vcf','-i',
    dest='input_vcf',
    help="input VCF file",
    type=str
  )
  parser.add_argument(
    '--idir',
    dest='idir',
    help="annotation summaries results dir",
    type=str
  )
  cargs =   parser.parse_args()


  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)
  VCF_NAME = Path(cargs.input_vcf).stem
  WDIR = os.getcwd()
  IDIR = cargs.idir



  print("Creating set list file")
  print("="*20)
  print(f"Config file: {os.path.basename(cargs.cfile)}")
  print(f"input VCF: {os.path.basename(cargs.input_vcf)}")
  print("="*20)





  main()