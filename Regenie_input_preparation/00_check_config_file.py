"""
Check the variables specified in the config file to make sure it has what it needs to run the scripts
"""
import yaml,os,shutil,sys,argparse
from collections import namedtuple
from itertools import compress



#%% General
# check all the programs exists
def __check_path_exists():
  """check to make sure the required paths exists
  """
  path_existence = {
    'bcftools': os.path.isfile(CONFIG.bcftools), 
    "tabix": os.path.isfile(CONFIG.tabix), 
    "apptainer": os.path.isfile(CONFIG.apptainer),
    "plink" : os.path.isfile(CONFIG.plink),
    "plink2" : os.path.isfile(CONFIG.plink2),
    "vep_docker": os.path.isfile(CONFIG.vep_docker_image),
    "script_dir": os.path.isdir(CONFIG.script_dir),
    "wdir": os.path.isdir(CONFIG.wdir),
    "helper_dir": os.path.isdir(os.path.join(CONFIG.script_dir,"python_helpers"))
  }
  assert(
    all(list(path_existence.values()))
  ),f"These are not found: {';'.join(list(compress(list(path_existence.keys()),[not x for x in list(path_existence.values())])))}"

  return


def __check_build():
  """Check to make sure the build is one that we work with (only GRCh38)
  ...might have to add other things to verify build (like check positions and what not...)
  """
  assert(
    CONFIG.genome_build == 'GRCh38'
  ),f"specified {CONFIG.genome_build}, not supported"
  return

def __check_names_def_consistencies():
  """checks the that all the categories for mask_names are defined in mask_definitions"""

  assert(
    len(list(CONFIG.mask_names.keys())) == 
    len(set(list(CONFIG.mask_names.keys())))
  ),f"duplicated studies in mask_names of config file"
  
  assert(
    len(list(CONFIG.mask_definitions.keys())) == 
    len(set(list(CONFIG.mask_definitions.keys())))
  ),f"duplicated studies in mask_definitions of config file"


  # all the mask names
  all_study_mask_names = {}
  for study,mask_names in CONFIG.mask_names.items():
    masks_tmp_list = []
    for i in list(mask_names.values()):
      masks_tmp_list += i.split(",")
      masks_tmp_list = list(set(masks_tmp_list))
    all_study_mask_names[study] = masks_tmp_list

  # all mask _definitions
  all_study_mask_def = {}
  for study,mask_def in CONFIG.mask_definitions.items():
    masks_tmp_list = []
    for i in list(mask_def.keys()):
      masks_tmp_list += i.split(",")
      masks_tmp_list = list(set(masks_tmp_list))
    all_study_mask_def[study] = masks_tmp_list
  
  assert(
    set(all_study_mask_def.keys()) == set(all_study_mask_names.keys())
  ),f"the studies are different between mask_definitions and mask_names"

  difference_checks = {}
  for study in all_study_mask_names.keys():
    difference_checks[study] = set(all_study_mask_names[study]) == set(all_study_mask_def[study])
  assert(
    all(list(difference_checks.values()))
  ),f"mask names and definitions are different for studies: {list(compress(list(difference_checks.keys()),[not x for x in difference_checks.values()]))}"

  return

def  __check_numeric_constant():
  """make sure the numeric constant values are either higher or lower
  """
  for k,v in CONFIG.CONST_NUMERIC.items():
    assert(
      v in ['higher','lower']
    ),f"Not 'higher' or 'lower' for Numeric constant plugin {k}"
  return

def __check_all_plugins_have_orders():
  """Check all plugins have specified orders in config file as required for aggregating annotations
  """
  all_plugins = set()
  for mask_def in CONFIG.mask_definitions.values():
    for plugins in mask_def.values():
      all_plugins.add(list(plugins.keys())[0])
  check_results = [a or b for a,b in zip([x in CONFIG.CONST for x in all_plugins],[ x in CONFIG.CONST_NUMERIC for x in all_plugins])]
  check_keys = all_plugins
  no_order_dict = dict(
    zip(check_keys,check_results)
  )
  assert(
    all(
      list(no_order_dict.values())
    )
  ),f"Plugins do not have order specified in config: {list(compress(list(no_order_dict.keys()),[not x for x in list(no_order_dict.values())]))}"

  return

def __check_def_have_plugin_orders():
  """each study defined plugins for generating annotataions. They should all have orders
  """
  for study,mask_def in CONFIG.mask_definitions.items():
    study_plugins = set(mask_def.keys())
    study_plugin_orders = CONFIG.plugin_orders[study]
    assert(
      study_plugins == set(study_plugin_orders)
    ),f"mask_definitions plugins and plugin orders are different for {study}"
  return

#%%
def main():
  __check_path_exists()
  __check_build()
  __check_names_def_consistencies()
  __check_numeric_constant()
  __check_all_plugins_have_orders()
  __check_def_have_plugin_orders()

  return

if __name__ == "__main__":
  parse = argparse.ArgumentParser()
  parse.add_argument(
    '--config_file','-c',
    dest='cfile',
    nargs=1,
    type=str,
    default = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"
  )
  cargs = parse.parse_args()

  # import mock
  # cargs = mock.Mock()
  # cargs.cfile = "/home/richards/kevin.liang2/scratch/exwas_pipeline/config/proj_config.yml"

  with open(cargs.cfile,'r') as ptr:
    params = yaml.safe_load(ptr)['proj_config']
  
  CONFIG = namedtuple('params',params.keys())(**params)


  main()