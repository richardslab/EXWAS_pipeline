def vcf_file = new File(params.annotation_vcf)
def vcf_file_name = vcf_file.getName()
process run_regenie_s2{
   publishDir("${params.outdir}/REGENIE_OUTPUTS/Regenie_S2",mode:"link")
  errorStrategy "terminate"
  label "big_memory","long"

  input:
    tuple val(regenie_input),val(ofile_suffix)
    val config_file
    val regenie_s1
    val nxtflow_genetic
    val nxtflow_genetic_type
    val nxtflow_annotation

  
  output:
    path "8_regenie_s2_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  set -o pipefail
  regenie_s2.py -c ${config_file} -i ${regenie_input} --regenie_s1 ${regenie_s1} --idir ${params.outdir}/REGENIE_INPUTS --nxtflow_genetic "${nxtflow_genetic}" --nxtflow_genetic_type ${nxtflow_genetic_type} --nxtflow_annotation "${nxtflow_annotation}" | tee 8_regenie_s2_${ofile_suffix}.logs
  """
}