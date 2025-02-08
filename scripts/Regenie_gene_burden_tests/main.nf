// This is the base file name produced by the scripts
// required so nextflow can check all outputs are produced
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

process run_regenie_s2{
  storeDir params.outdir

  input:
    tuple val(regenie_input),val(ofile_suffix)
    val config_file
    val wdir
    val nxtflow_genetic
    val nxtflow_genetic_type
    val nxtflow_annotation
    path "7_regenie_s1.logs"

  
  output:
    path "8_regenie_s2_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  set -o pipefail

  python -u ${baseDir}/modules/Regenie_gene_burden_tests/02_regenie_s2.py -c ${config_file} -i ${regenie_input} --wdir ${wdir} --nxtflow_genetic "${nxtflow_genetic}" --nxtflow_genetic_type ${nxtflow_genetic_type} --nxtflow_annotation "${nxtflow_annotation}" | tee 8_regenie_s2_${ofile_suffix}.logs
  """
}