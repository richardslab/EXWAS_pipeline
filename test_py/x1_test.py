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
  parser = argparse.ArgumentParser()
  parser.add_argument('-c',type=str,nargs=1,dest='cfile')
  cargs = parser.parse_args()
  print(cargs.cfile)
  return

def print_env():
  sp.run(
    ['which','python'],check=True
  )
  return
print_env()

