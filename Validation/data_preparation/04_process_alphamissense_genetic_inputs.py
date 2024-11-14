import os,shutil,yaml
import subprocess as sp

pgen_file = "/scratch/richards/yiheng.chen/UKB_data/WES_UKB_data/ukb_merged_1-22"
OUTDIR = "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/alphamissense_results"
ofile = "ukb_merged_1-22_sitesonly.vcf"
PLINK="/home/richards/kevin.liang2/scratch/programs/plink2"
plink_cmd = [
  PLINK,
  '--pfile',pgen_file,
  '--make-just-pvar',"cols=xheader,vcfheader,qual,filter,info",
  '--output-chr','chrM',
  '--out',
  os.path.join(OUTDIR,ofile)
]
res = sp.run(plink_cmd,check=True)
assert(res.returncode == 0)