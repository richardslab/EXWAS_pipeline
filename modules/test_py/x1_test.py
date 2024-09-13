import subprocess as sp
import argparse

def hello():
  print("hi")
  return

def write_hello():
  sp.run(
    ['touch','/Users/kevinliang/Desktop/work/working/exwas_pipelines/results/test_output/hi2.txt'],check=True
  )
  return

def get_param():
  # sp.run(
  #   'echo $CONDA_PREFIX',check=True,shell=True
  # )
  sp.run(
    'echo $(which python)',check=True,shell=True
  )
  # sp.run(
  #   'conda info',check=True,shell=True
  # )
  # sp.run(
  #   'echo $PATH',check=True,shell=True
  # )
  import pyreadr
  parser = argparse.ArgumentParser()
  parser.add_argument('-c',type=str,nargs=1,dest='cfile')
  cargs = parser.parse_args()
  print(cargs.cfile)
  return

get_param()

