def vcf_file = new File(params.annotation_vcf)
def vcf_file_name = vcf_file.getName()
process run_regenie_s1{
   storeDir params.outdir

  input:
    val config_file
    val wdir
  
  output:
    path "Regenie_S1/7_Regenie_S1_pred.list"
    path "7_regenie_s1.logs",emit: "log"
  
  script:
  """
  set -o pipefail

  python -u ${baseDir}/modules/Regenie_gene_burden_tests/01_regenie_s1.py -c ${config_file} --wdir ${wdir} | tee 7_regenie_s1.logs
  """
}