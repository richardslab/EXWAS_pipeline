// This is the base file name produced by the scripts
// required so nextflow can check all outputs are produced
def vcf_file = new File(params.annotation_vcf)
def vcf_file_name = vcf_file.getName()

process run_regenie_s1{
   storeDir params.outdir

  input:
    val config_file
    val wdir
    tuple val(regenie_input) val(ofile_suffix)
  
  output:
    path "1_regenie_S1_OUT_pred_${ofile_suffix}.list"
    path "1_regenie_s1_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  python -u ${baseDir}/modules/Regenie_gene_burden_tests/01_regenie_s1.py -c ${config_file} -i ${regenie_input} --wdir ${wdir} | tee 1_regenie_s1_${ofile_suffix}.logs
  """
}

process run_regenie_s2{
   storeDir params.outdir

  input:
    val config_file
    val wdir
    tuple val(regenie_input) val(ofile_suffix)
    val file_type
    path "1_regenie_s1_${ofile_suffix}.logs"
  
  output:
    path "2_regenie_s2_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  python -u ${baseDir}/modules/Regenie_gene_burden_tests/02_regenie_s2.py -c ${config_file} -i ${regenie_input} --wdir ${wdir} --regenie_file_type ${file_type} | tee 2_regenie_s2_${regenie_input}.logs
  """
}