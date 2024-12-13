""" Create the annotation file for each study based on the masks defined in the config file.

# parse the VEP anotations so we have all variant annotations for each gene with consequence sorted.
# put in sqlite3 database and index it for faster access
# would be similar to using multiple files storing all variants per annotation
"""

import shutil,os,sys,yaml,pyreadr,re,json,gzip,sqlite3
import pandas as pd
import numpy as np
import subprocess as sp
from collections import namedtuple
import argparse
from pathlib import Path

def __init_annotation_db(db_file):
  try:
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    cur.execute(
    """
      CREATE TABLE vep_summaries(SNP,location,plugin,plugin_consequence,gene,var_consequence)
    """
    )
  except Exception as e:
    print(f"SQLITE3 error: {e}")
    conn.commit()
    conn.close()
    raise
  return

def __obtain_annotation_file_headers(expected_annotation_file):
  """Obtain the header of the file.
  Do everything separately so can parse the annotation file in parallele and not line by line
  Args:
      expected_annotation_file (str): the path to annotation file
  """
  n_header_rows = 0
  with gzip.open(expected_annotation_file,'rt') as ptr:
    headers = []
    id_idx = None
    location_idx = None
    gene_idx = None
    annotation_idx = None
    for line in ptr:
      n_header_rows += 1
      if re.match("## ",line):
        continue
      # this is the header line
      if re.match("#Uploaded_variation",line):
        headers = [x.strip() for x in line.split()]
        assert(
          CONFIG.vep_variant_ID in headers and 
          CONFIG.vep_variant_location in headers and 
          CONFIG.vep_gene in headers and 
          CONFIG.vep_annotations in headers and
          CONFIG.vep_var_consequence in headers
        )
        id_idx = headers.index(CONFIG.vep_variant_ID)
        location_idx = headers.index(CONFIG.vep_variant_location)
        gene_idx = headers.index(CONFIG.vep_gene)
        annotation_idx = headers.index(CONFIG.vep_annotations)
        vep_var_consequence_idx = headers.index(CONFIG.vep_var_consequence)
        continue
      # finished all the headers, didn't find the header line, then there is a problem
      # with the annotation file
      if not re.match("#",line):
        if (id_idx == None or location_idx == None or gene_idx == None or annotation_idx == None or vep_var_consequence_idx == None or len(headers) == 0):
          assert(False),"did not find header or relevant columns in {expected_annotation_file}"
        else:
          break
  return n_header_rows,headers,id_idx,location_idx,gene_idx,annotation_idx,vep_var_consequence_idx

def __get_consequence_summaries(line,id_idx,gene_idx,location_idx,annotation_idx,vep_var_consequence_idx):
  line_elem = [x.strip() for x in line.split("\t")]
  variant_id = line_elem[id_idx]
  gene = line_elem[gene_idx]
  vep_var_consequence = line_elem[vep_var_consequence_idx]
  if gene == '-':
    return None,None
  location = line_elem[location_idx]
  annotation = line_elem[annotation_idx]  
  var_plugin_consequence_summaries = parse_vep.parse_var_consequence(annotation,ALL_CONST)
  var_info = {"variant_id" : variant_id,"gene":gene,"location":location,"var_consequence":vep_var_consequence}
  return var_plugin_consequence_summaries,var_info

def __get_update_info(cur,plugin,new_consequence,var_info):
  # based on query
  query_location_idx = 0
  query_consequence_idx = 1
  update_info = {
    "SNP":var_info['variant_id'],
    "location":var_info['location'],
    "plugin":None,
    "plugin_consequence":None,
    "gene":var_info['gene'],
    "var_consequence":var_info['var_consequence'],
    "update":False
  }
  query_params={k:var_info[k] for k in ['gene','variant_id']}
  query_params.update({"plugin":plugin})
  gene_var_query = cur.execute(
  """
    SELECT location,plugin_consequence from vep_summaries WHERE gene == :gene AND SNP == :variant_id AND plugin == :plugin
  """,query_params
  )
  gene_vars = gene_var_query.fetchall()
  # for each plugin, each variant, each gene, there should be at most 1 entry
  assert(len(gene_vars) <= 1),f"duplicated variant for same gene {var_info['gene']},{var_info['variant_id']},{plugin}"
  # update consequence if more severe
  if len(gene_vars) == 1:
    update_info['plugin'] = plugin
    # check location is the same for same variant
    gene_vars = gene_vars[0] # gene vars initially a tuple
    assert(
      gene_vars[query_location_idx] == var_info['location']
    ),f"inconsistent annotation var: {var_info['variant_id']}. Found {gene_vars[query_location_idx]}. {query_location_idx} {gene_vars}"
    # update consequence to keep most severe across different transcripts
    if new_consequence == None:
      return update_info
    if gene_vars[query_consequence_idx] == None:
      update_info['plugin'] = plugin
      update_info['plugin_consequence'] = new_consequence
      update_info['update'] =True
      return update_info
    if plugin in CONSTANT:
      existing_consequence_priority = CONSTANT[plugin].index(gene_vars[query_consequence_idx])
      new_consequence_priority = CONSTANT[plugin].index(new_consequence)
      if new_consequence_priority < existing_consequence_priority:
        update_info['plugin'] = plugin
        update_info['plugin_consequence'] = new_consequence
        update_info['update'] =True
        return update_info
    elif plugin in CONST_NUMERIC:
      if CONST_NUMERIC[plugin] == 'higher':
        if new_consequence > gene_vars[query_consequence_idx]:
          update_info['plugin'] = plugin
          update_info['plugin_consequence'] = new_consequence
          update_info['update'] =True
          return update_info
      elif CONST_NUMERIC[plugin] == 'lower':
        if new_consequence < gene_vars[query_consequence_idx]:
          update_info['plugin'] = plugin
          update_info['plugin_consequence'] = new_consequence
          update_info['update'] =True
          return update_info
      else:
        assert(False),f"unknown numeric plugin orders {plugin}"
    else:
      assert(False),f"Unknown plugin type {plugin}"
  else:
    update_info['plugin'] = plugin
    update_info['plugin_consequence'] = new_consequence
  return update_info

def __update_annotation_db(cur,update_info):
  # constant for indexing database results
  fields_to_update = {k:v for k,v in update_info.items() if v != None}
  # if all fields are new == new entry
  if not fields_to_update['update'] and 'plugin_consequence' in fields_to_update:
    cur.execute(
      """
      INSERT INTO vep_summaries VALUES(:SNP,:location,:plugin,:plugin_consequence,:gene,:var_consequence) 
      """,fields_to_update
    )
  elif fields_to_update['update']:
    assert('plugin_consequence' in fields_to_update),"what are we updating..."
    cur.execute(
      """
      UPDATE vep_summaries SET plugin_consequence = :plugin_consequence, var_consequence = :var_consequence WHERE SNP = :SNP AND plugin = :plugin AND gene = :gene AND location = :location
      """,fields_to_update
    )
  return

def __index_db(db_file,reindex=False):
  try:  
    # index the table
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    cur.execute(
    """
    CREATE UNIQUE INDEX gene_idx ON vep_summaries(gene,SNP,location,plugin);
    """
    )
    cur.execute(
    """
    CREATE UNIQUE INDEX plugin_var_consequence ON vep_summaries(plugin,plugin_consequence,var_consequence,SNP,gene)
    """
    )
  except sqlite3.OperationalError as e:
    print(f"SQLITE3 error: {e}")
    conn.close()
    raise
  return

def main():
  expected_annotation_file = os.path.join(WDIR,f'3_annotation_results_{VCF_NAME}.txt.gz')
  db_file = os.path.join(
    WDIR,f"5_1_vep_summaries_{VCF_NAME}.sqlite3.db"
  )
  
  print("Summarizing all VEP annotations")
  print("*"*20)
  print("The most severe consequence per gene (across transcripts) are kept")
  print(f"Consequence order obtained from {cargs.cfile}")
  print(f"Summaries are stored in {db_file}")
  print("")

  # initialize database and create index
  __init_annotation_db(db_file)
  __index_db(db_file)

  # obtain headers
  n_header_rows,headers,id_idx,location_idx,gene_idx,annotation_idx,vep_var_consequence_idx = __obtain_annotation_file_headers(expected_annotation_file)
  gunzip_cmd = sp.Popen(
    [
      'gunzip','-c',expected_annotation_file
    ],stdout=sp.PIPE
  )
  n_rows = sp.run(
    ['wc','-l'],
    stdin = gunzip_cmd.stdout,check=True,capture_output=True
  ).stdout.decode('utf-8').split()[0]
  gunzip_cmd.stdout.close()
  n_rows = float(n_rows)
  print(f"{expected_annotation_file} have {n_rows} variants")
  
  # update the files line by line and writes to db in chunks
  # writes to db every <chunk_size> lines
  chunk_size = 5000
  try:
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    with gzip.open(expected_annotation_file,'rt') as ptr:
      line_ct = 0
      # whether to commit the db (i.e., write to file)
      write_chunks = False
      for line in ptr:
        line_ct += 1
        if line_ct % 1000 == 0:
          print(f"{line_ct}/{n_rows}: {'{:.3e}'.format(line_ct/n_rows*100)}% completed")
        if line_ct % chunk_size == 0:
          write_chunks = True
        if bool(re.match("#",line)):
          continue
        var_plugin_consequence_summaries,var_info = __get_consequence_summaries(line,id_idx,gene_idx,location_idx,annotation_idx,vep_var_consequence_idx)
        if var_plugin_consequence_summaries is None:
          continue
        for plugin,new_consequence in var_plugin_consequence_summaries.items():
          update_info = __get_update_info(cur,plugin,new_consequence,var_info)
          __update_annotation_db(cur,update_info)
        # write and index after each
        if write_chunks:
          conn.commit()
          write_chunks = False
    conn.close()
  except Exception as e:
    print(f"SQLITE3 error: {e}")
    conn.close()
    raise

  try:
    conn = sqlite3.connect(db_file)
    cur = conn.cursor()
    # it is a tuple
    distinct_genes = cur.execute("SELECT count(*) FROM (SELECT distinct gene FROM vep_summaries)").fetchone()
    distinct_genes = distinct_genes[0]
    distinct_variants = cur.execute("SELECT count(*) FROM (SELECT distinct SNP FROM vep_summaries)").fetchone()
    distinct_variants = distinct_variants[0]
  except Exception as e:
    print(f"SQLITE3 ERROR: {3}")
    con.close()
    res = None
    raise

  print(f"Parsed {distinct_variants} variants")
  print(f"generated annotations for {distinct_genes} genes")
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
    cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/alphamis_exact_configs/proj_config.yml"
    cargs.wdir="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/alphamiss_exact"
    cargs.input_vcf="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/alphamissense_results/ukb_merged_1-22_sitesonly.vcf"
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
  ALL_CONST = {}
  ALL_CONST.update(CONSTANT)
  ALL_CONST.update(CONST_NUMERIC)
  
  sys.path.append(os.path.dirname(__file__))
  from python_helpers.vep_helpers import parse_vep
  from python_helpers.vep_helpers import parse_vep_headers



  main()