""" Create the annotation file for each study based on the masks defined in the config file.

# parse the VEP anotations so we have all variant annotations for each gene with consequence sorted.
# put in sqlite3 database and index it for faster access
# would be similar to using multiple files storing all variants per annotation
"""

import shutil,os,sys,yaml,pyreadr,re,json,gzip,sqlite3
import pandas as pd
import numpy as np
from collections import namedtuple
import argparse
from pathlib import Path

def __init_annotation_db(db_file):
  conn = sqlite3.connect(db_file)
  cur = conn.cursor()
  cur.execute(
  """
    CREATE TABLE vep_summaries(SNP,location,plugin,consequence,gene)
  """
  )
  conn.commit()
  conn.close()
  return

def __update_annotation_db(db_file,variant_id,gene,location,
var_consequence_summaries):
  # constant for indexing database results
  location_idx = 0
  consequence_idx = 1
  conn = sqlite3.connect(db_file)
  cur = conn.cursor()
  for plugin,new_consequence in var_consequence_summaries.items():
    update_info = {
      "SNP":None,
      "location":None,
      "plugin":None,
      "consequence":None,
      "gene":None
    }
    gene_var_query = cur.execute(
    """
      SELECT location,consequence from vep_summaries WHERE gene == :gene AND SNP == :variant_id AND plugin == :plugin
    """,{"gene":gene,"variant_id":variant_id,"plugin":plugin}
    )
    gene_vars = gene_var_query.fetchall()
    # for each plugin, each variant, each gene, there should be at most 1 entry
    assert(len(gene_vars) <= 1),f"duplicated variant for same gene {gene},{variant_id},{plugin}"
    # update consequence if more severe
    if len(gene_vars) == 1:
      update_info['plugin'] = plugin
      # check location is the same for same variant
      gene_vars = gene_vars[0]
      assert(
        gene_vars[location_idx] == location
      ),f"inconsistent annotation var: {variant_id}"
      # update consequence to keep most severe across different transcripts
      if new_consequence == None:
        continue
      if gene_vars[consequence_idx] == None:
        update_info['new_consequence'] = new_consequence
        continue
      if plugin in CONSTANT:
        existing_consequence_priority = CONSTANT[plugin].index(gene_vars[consequence_idx])
        new_consequence_priority = CONSTANT[plugin].index(new_consequence)
        if new_consequence_priority < existing_consequence_priority:
          update_info['new_consequence'] = new_consequence
          continue
      elif plugin in CONST_NUMERIC:
        if CONST_NUMERIC[plugin] == 'higher':
          if new_consequence > gene_vars[consequence_idx]:
            update_info['new_consequence'] = new_consequence
            continue
        elif CONST_NUMERIC[plugin] == 'lower':
          if new_consequence < gene_vars[consequence_idx]:
            update_info['new_consequence'] = new_consequence
            continue
        else:
          assert(False),f"unknown numeric plugin orders {plugin}"
      else:
        assert(False),f"Unknown plugin type {plugin}"
    else:
      update_info['SNP'] = variant_id
      update_info['location'] = location
      update_info['plugin'] = plugin
      update_info['consequence'] = new_consequence
      update_info['gene'] = gene
    
    fields_to_update = {k:v for k,v in update_info.items() if v != None}
    # if all fields are new == new entry
    if len(fields_to_update) == 5:
      cur.execute(
        """
        INSERT INTO vep_summaries VALUES(:SNP,:location,:plugin,:consequence,:gene) 
        """,fields_to_update
      )
    elif "consequence" in fields_to_update:
      cur.execute(
        """
        UPDATE vep_summaries SET consequence = :consequence WHERE SNP = :SNP AND plugin = :plugin AND gene = :gene AND location = :location
        """,fields_to_update
      )
  conn.commit()
  conn.close()

  return

def main():
  expected_annotation_file = os.path.join(WDIR,f'3_annotation_results_{VCF_NAME}.txt')
  db_file = os.path.join(
    WDIR,f"5_1_vep_summaries_{VCF_NAME}.sqlite3.db"
  )
  
  print("Summarizing all VEP annotations")
  print("*"*20)
  print("The most severe consequence per gene (across transcripts) are kept")
  print(f"Consequence order obtained from {cargs.cfile}")
  print(f"Summaries are stored in {db_file}")
  __init_annotation_db(db_file)
  print("")

  line_num = 0
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
      # this is the header line
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
    
      # finished all the headers, didn't find the header line, then there is a problem
      # with the annotation file
      if not re.match("#",line) and (id_idx == None or location_idx == None or gene_idx == None or annotation_idx == None or len(headers) == 0):
        assert(False),"did not find header or relevant columns in {expected_annotation_file}"

      line_elem = [x.strip() for x in line.split()]
      variant_id = line_elem[id_idx]
      gene = line_elem[gene_idx]
      location = line_elem[location_idx]
      annotation = line_elem[annotation_idx]  
      var_consequence_summaries = parse_vep.parse_var_consequence(annotation,CONSTANT)
      __update_annotation_db(db_file,variant_id,gene,location,var_consequence_summaries)
  
  # index the table
  conn = sqlite3.connect(db_file)
  cur = conn.cursor()
  cur.execute(
  """
  CREATE UNIQUE INDEX gene_index ON vep_summaries(gene,SNP,plugin,consequence)
  """
  )
  cur.execute(
  """
  CREATE UNIQUE INDEX plugin_index ON vep_summaries(plugin,consequence,SNP,gene)
  """
  )
  cur.execute(
  """
  CREATE UNIQUE INDEX var_index ON vep_summaries(SNP,gene,plugin,consequence)
  """
  )
  res = cur.execute("SELECT count(*) FROM (SELECT distinct gene FROM vep_summaries)").fetchone()
  res = res[0]
  conn.commit()
  conn.close()
  print(f"generated annotations for {res} genes")
  print("="*20)
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
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/sitesonly_VCF/wes_qc_chr10_sitesonly.vcf"
    __file__ = "/home/richards/kevin.liang2/scratch/exwas_pipeline/src/modules/Regenie_input_preparation/04_1_create_annotation_summaries.py"
    print("TEST")


  assert(os.path.isfile(cargs.cfile)),'config file is missing'
  assert(os.path.isfile(cargs.input_vcf)),'input vcf is missing'
  assert(cargs.wdir),'output directory missing'
  print("Creating annotation summaries")
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
  CONSTANT = CONFIG.CONST
  CONST_NUMERIC = CONFIG.CONST_NUMERIC
  
  sys.path.append(os.path.dirname(__file__))
  from python_helpers.vep_helpers import parse_vep
  from python_helpers.vep_helpers import parse_vep_headers



  main()