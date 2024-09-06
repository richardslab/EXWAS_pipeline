""" Create the annotation file for each study based on the masks defined in the config file.

# parse the VEP anotations so we have all variant annotations for each gene with consequence sorted.
"""

import shutil,os,sys,yaml,pyreadr,re
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse

def main():
  expected_annotation_file = os.path.join(CONFIG.wdir,f'2_{VCF_NAME}_vcf_final_annotation.txt')

  outdir = os.path.join(CONFIG.wdir,'vep_consequence_summaries')
  os.makedirs(outdir,exist_ok=True)
  
  vep_summaries = dict()
  with open(expected_annotation_file,'r') as ptr:
    headers = []
    id_idx = None
    location_idx = None
    gene_idx = None
    annotation_idx = None
    for line in ptr:
      if re.match("## ",line):
        continue
      if re.match("#Uploaded_variation",line):
        headers = [x.strip() for x in line.split()]
        assert(
          all(
            CONFIG.vep_variant_ID in headers and 
            CONFIG.vep_variant_location in headers and 
            CONFIG.vep_gene in headers and 
            CONFIG.vep_annotations in headers 
          )
        )
        id_idx = headers.index(CONFIG.vep_variant_ID)
        location_idx = headers.index(CONFIG.vep_variant_location)
        gene_idx = headers.index(CONFIG.vep_gene)
        annotation_idx = headers.index(CONFIG.vep_annotations)
        continue
      # finished all the headers, didn't find what we needed
      if not re.match("#",line) and (id_idx == None or location_idx == None or gene_idx == None or annotation_idx == None):
        assert(False),"did not find header or relevant columns"

      if len(headers > 0):
        line_elem = [x.strip() for x in line.split()]
        variant_id = line_elem[id_idx]
        gene = line_elem[gene_idx]
        location = line_elem[location_idx]
        annotation = line_elem[annotation_idx]  
        var_consequence_summaries = parse_vep.parse_var_consequence(annotation)
        if gene in vep_summaries:
          vep_summaries[gene]+={
            "variant_id":variant_id,
            "location":location,
            "consequences":var_consequence_summaries
          }
      

  return

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '--config_file','-c',
    dest='cfile',
    default="/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml",
    help='configuration yaml file'
  )
  cargs =   parser.parse_args()

  # import mock
  # cargs = mock.Mock()
  # cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  print(f"Using {os.path.basename(cargs.cfile)}")



  with open(cargs.cfile,'r') as ptr:
    params = yaml.full_load(ptr)['proj_config']
  CONFIG = namedtuple("params",params.keys())(**params)

  VCF_NAME = os.path.basename(CONFIG.input_vcf)
  sys.path.append(CONFIG.script_dir)
  from python_helpers.vep_helpers import parse_vep


  main()