""" Create the annotation file for each study based on the masks defined in the config file.

# parse the VEP anotations so we have all variant annotations for each gene with consequence sorted.
"""

import shutil,os,sys,yaml,pyreadr,re,json,gzip
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse

def main():
  expected_annotation_file = os.path.join(WDIR,f'2_{VCF_NAME}_vcf_final_annotation.txt')

  outdir = os.path.join(WDIR,'vep_consequence_summaries')
  os.makedirs(outdir,exist_ok=True)
  
  line_num = 0
  vep_summaries = dict()
  with open(expected_annotation_file,'r') as ptr:
    headers = []
    id_idx = None
    location_idx = None
    gene_idx = None
    annotation_idx = None
    for line in ptr:
      line_num += 1
      if re.match("## ",line):
        continue
      if re.match("#Uploaded_variation",line):
        headers = [x.strip() for x in line.split()]
        assert(
          CONFIG.vep_variant_ID in headers and 
          CONFIG.vep_variant_location in headers and 
          CONFIG.vep_gene in headers and 
          CONFIG.vep_annotations in headers 
        )
        id_idx = headers.index(CONFIG.vep_variant_ID)
        location_idx = headers.index(CONFIG.vep_variant_location)
        gene_idx = headers.index(CONFIG.vep_gene)
        annotation_idx = headers.index(CONFIG.vep_annotations)
        continue
    
      # finished all the headers, didn't find what we needed
      if not re.match("#",line) and (id_idx == None or location_idx == None or gene_idx == None or annotation_idx == None):
        assert(False),"did not find header or relevant columns"

      if len(headers) > 0:
        line_elem = [x.strip() for x in line.split()]
        variant_id = line_elem[id_idx]
        gene = line_elem[gene_idx]
        location = line_elem[location_idx]
        annotation = line_elem[annotation_idx]  
        var_consequence_summaries = parse_vep.parse_var_consequence(annotation,CONSTANT)
        if gene in vep_summaries:
          # if same variant, update consequence if it is more severe
          if variant_id in vep_summaries[gene]:
            assert(
              location == vep_summaries[gene][variant_id]['location']
            ),f"inconsistent annotation var: {variant_id} gene: {gene} line {line_num}"
            existing_consequence = vep_summaries[gene][variant_id]['consequences']
            for k,v in var_consequence_summaries.items():
              new_value = var_consequence_summaries[k]
              
              # if this plugin was not used before, then replace with new value
              # if this plugin is not used now, then move on
              if k not in existing_consequence:
                vep_summaries[gene][variant_id]['consequences'][k]=new_value
                continue
              
              old_value = existing_consequence[k]
              if new_value == None:
                continue
              if old_value == None:
                if new_value != None:
                  vep_summaries[gene][variant_id]['consequences'][k]=new_value
                continue
              # for categorical plugins (e.g., LoF, IMPACT, etc)
              if k in CONSTANT:
                if k in existing_consequence:
                  old_consequence_priority = CONSTANT[k].index(old_value)
                  new_consequence_priority = CONSTANT[k].index(new_value)
                  if new_consequence_priority < old_consequence_priority:
                    vep_summaries[gene][variant_id]['consequences'][k] = v
              elif k in CONST_NUMERIC:
                if k in existing_consequence:
                  if CONST_NUMERIC[k] == 'higher':
                    if new_value > old_value:
                      vep_summaries[gene][variant_id]['consequences'][k] = v
                  elif CONST_NUMERIC[k] == 'lower':
                    if new_value < old_value:
                      vep_summaries[gene][variant_id]['consequences'][k] = v
                  else:
                    assert(False),f"unknown numeric interpretation {k}"
              else:
                assert(False),f"unknown plugin {k}"
          else:
            vep_summaries[gene] = {
              variant_id:{
                "location":location,
                "consequences":var_consequence_summaries
              }
            }
        else:
          vep_summaries[gene] = {
            variant_id:{
              "location":location,
              "consequences":var_consequence_summaries
            }
          }
      else:
        assert(False),f"should not reach here!"

  vep_summarie_file = os.path.join(
    outdir,f"{VCF_NAME}_vep_summaries.json.gz"
  )
  with gzip.open(vep_summarie_file,'wt',encoding='utf-8') as ptr:
    json.dump(
      vep_summaries,
      ptr,ensure_ascii=False
    )
  print(f"generated annoations for {len(vep_summaries.keys())} genes")
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
  CONSTANT = CONFIG.CONST
  CONST_NUMERIC = CONFIG.CONST_NUMERIC
  
  sys.path.append(CONFIG.script_dir)
  from python_helpers.vep_helpers import parse_vep
  from python_helpers.vep_helpers import parse_vep_headers



  main()