"""
  Add 'chr' prefix to variant names
"""
import os,shutil,re,glob
import subprocess as sp

# constants
PLINK="./plink2"
OUTDIR="./exwas_data"
INDIR="./wes_qc_chr"

chrom_idx = 0
pos_idx = 1
id_idx = 2
ref_idx = 3
alt_idx = 4
def fix_chr_notations(all_pgen_files):
  for each_f in all_pgen_files:
    n_line_original = float(sp.run(
      ['wc','-l',f"{each_f}.pvar"],check=True,capture_output=True
    ).stdout.decode('utf-8').split()[0])
    ofile = os.path.join(OUTDIR,f"{os.path.basename(each_f)}")
    with open(f"{each_f}.pvar",'r') as iptr, open(f"{ofile}.pvar",'w') as optr:
      for line in iptr:
        line = line.strip()
        if re.match("#",line):
          res = optr.write(f"{line}\n")
        else:
          line_elem = line.split("\t")
          assert(len(line_elem) == 5)
          n_chr = f"chr{line_elem[chrom_idx]}"
          n_id = f"chr{line_elem[id_idx]}"
          res = optr.write("\t".join([n_chr,line_elem[pos_idx],n_id,line_elem[ref_idx],line_elem[alt_idx]]) + "\n")
    n_line_new = float(sp.run(
    ['wc','-l',f"{ofile}.pvar"],check=True,capture_output=True
  ).stdout.decode('utf-8').split()[0])
  assert(n_line_new == n_line_original),f"{each_f}"
  return

def rename_file(all_sitesonly_vcf):
  for each_f in all_sitesonly_vcf:
    os.rename(each_f,re.sub("\\.pvar","",each_f))
  return

def main():
  all_pgen_files = [re.sub("\\.psam","",x) for x in glob.glob(os.path.join(INDIR,"*.psam"))]
  # convert to sites only VCF
  fix_chr_notations(all_pgen_files)
  
  
  return

if __name__ == "__main__":
  main()