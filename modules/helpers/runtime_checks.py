import os
def check_runtime_environment():
  """check runtime environment

  For some reason (have to double check with nextflow people), that running nextflow with conda environment enabled that activates a different conda environemnt causes issues (the path variable order is messed up)

  Running with no conda environment is fine, so just make sure that is the cases
  """
  conda_env = os.environ['CONDA_PREFIX']
  assert conda_env == '',f"Deactivate conda before running nextflow. Found {conda_env}"
  return
check_runtime_environment()