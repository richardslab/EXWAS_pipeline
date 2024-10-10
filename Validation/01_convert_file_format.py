"""
  Convert PGEN files to sites only VCF files
"""
import os,shutil,re,glob
import subprocess as sp

# constants
PLINK="/home/richards/kevin.liang2/scratch/programs/plink2"
OUTDIR="/home/richards/kevin.liang2/scratch/exwas_pipeline/results/sitesonly_VCF"
INDIR="/scratch/richards/guillaume.butler-laporte/ukb_covid/wes_qc_chr"

def convert_to_sites_only(all_pgen_files):
  for each_f in all_pgen_files:
    plink_cmd = [
      PLINK,
      '--pfile',each_f,
      '--make-just-pvar',"cols=xheader,vcfheader,qual,filter,info",
      '--out',
      os.path.join(OUTDIR,f"{os.path.basename(each_f)}_sitesonly.vcf")
    ]
    res = sp.run(plink_cmd,check=True)
    assert(res.returncode == 0)
  return

def rename_file(all_sitesonly_vcf):
  for each_f in all_sitesonly_vcf:
    os.rename(each_f,re.sub("\\.pvar","",each_f))
  return

def main():
  all_pgen_files = [re.sub("\\.psam","",x) for x in glob.glob(os.path.join(INDIR,"*.psam"))]
  # convert to sites only VCF
  convert_to_sites_only(all_pgen_files)
  # rename file extensions because it is --make-just-pvar so pvar output
  all_sitesonly_vcf = glob.glob(os.path.join(OUTDIR,"*_sitesonly.vcf.pvar"))
  rename_file(all_sitesonly_vcf)
  
  return

if __name__ == "__main__":
  main()