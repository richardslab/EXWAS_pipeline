""" Create the set list file, which is shared by all studies. It is based on the input VCF.

"""

import shutil,os,sys,yaml,pyreadr,re,json,gzip
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse
from pathlib import Path

def main():
  vep_summarie_file = os.path.join(
    WDIR,'vep_consequence_summaries',f"{VCF_NAME}_vep_summaries.json.gz"
  )
  assert(
    os.path.isfile(vep_summarie_file)
  ),f"missing summary file"

  with gzip.open(vep_summarie_file,'r') as ptr:
    vep_summaries = json.load(ptr)

  setlist_file = os.path.join(
    WDIR,'vep_consequence_summaries',f"{VCF_NAME}.setlist"
  )
  with open(setlist_file,'w') as ptr:
    for gene,vars in vep_summaries.items():
      min_position = None
      gene_chr = None
      all_vars = []
      for var,var_info in vars.items():
        all_vars.append(var)
        var_location = var_info['location']
        var_location_info = var_location.split(":")
        assert(
          len(var_location_info) == 2
        ),f"unknown var notation {var_location}"
        if (re.search("-",var_location_info[1])):
          var_pos = float(var_location_info[1].split("-")[0])
        else:  
          var_pos = float(var_location_info[1])
        var_chr = str(re.sub("chr","",var_location_info[0]))
        if min_position == None:
          min_position == var_pos
        else:
          if var_pos < min_position:
            min_position = var_pos
        if gene_chr == None:
          gene_chr = var_chr
        else:
          assert(gene_chr == var_chr),f"different chr for same gene {gene} {var}"
      all_vars_str = ",".join(all_vars)
      ptr.write(
        "{gene}\t{gene_chr}\t{min_position}\t{all_vars_str}\n"
      )


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
  sys.path.append(CONFIG.Regenie_input_prep_scripts)




  main()