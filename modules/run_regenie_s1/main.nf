def vcf_file = new File(params.annotation_vcf)
def vcf_file_name = vcf_file.getName()
process run_regenie_s1{
   publishDir("${params.outdir}/REGENIE_OUTPUTS/Regenie_S1",mode:"link")
  errorStrategy "terminate"
  label "big_memory","long"

  input:
    val config_file
  
  output:
    path "7_Regenie_S1_pred.list", emit: "s1_output"
    path "7_regenie_s1.logs",emit: "log"
  
  script:
  """
  set -o pipefail

  regenie_s1.py -c ${config_file} | tee 7_regenie_s1.logs
  """
}