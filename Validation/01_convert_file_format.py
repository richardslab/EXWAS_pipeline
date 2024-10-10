"""
  Convert PGEN files to sites only VCF files
"""
import os,shutil,re,glob
import subprocess as sp

# constants
PLINK=""
OUTDIR=""
INDIR=""

def main():
  
  all_pgen_files = [re.sub("\\.psam","",x) for x in glob.glob(os.path.join(INDIR,"*.psam"))]
  for each_f in all_pgen_files:
    plink_cmd = [
      PLINK,
      '--pgen',each_f,
      '--make-just-pvar','pvar-cols','xheader,vcfheader,qual,filter,info',
      '--export','bgz','vcf',
    ]
    res = sp.run(plink_cmd,check=True)
    assert(res.returncode == 0)

  return

if __name__ == "__main__":
  main()