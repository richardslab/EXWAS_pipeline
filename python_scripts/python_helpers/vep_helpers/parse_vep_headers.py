import re

def get_vep_plugins(vep_file):
  """
    Helper to get plugin names from VEP outputs

    # obtain the plugins used in vep
    ## the extra column key is the last set of info before the results
    ## https://useast.ensembl.org/info/docs/tools/vep/vep_formats.html
  """
  start_serach_str = "## Extra column keys:"
  end_search_str="#Uploaded_variation"

  extra_columns = []
  with open(vep_file,'r') as ptr:
    start_search = False
    end_search = False
    for line in ptr:
      if re.match(start_serach_str,line):
        start_search=True
        continue
      if re.match(end_search_str,line):
        end_search = True
        break
      if start_search and not end_search:
        extra_columns += [re.sub("##","",line.split(":")[0]).strip()]
  assert(len(extra_columns) > 0),f"No plugin found"
  return extra_columns