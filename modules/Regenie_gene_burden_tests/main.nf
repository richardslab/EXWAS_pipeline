// This is the base file name produced by the scripts
// required so nextflow can check all outputs are produced
def vcf_file = new File(params.input_vcf)
def vcf_file_name = vcf_file.getName()

process run_regenie_s1{
   storeDir params.outdir

  input:
    val config_file
    val input_vcf
    val wdir
    path "6_create_setlist_file.logs"
  
  output:
    path "7_${vcf_file_name}_regenie_S1_OUT.pred"
    path "7_regenie_s1.logs",emit: "log"
  
  script:
  """
  python -u ${baseDir}/modules/Regenie_gene_burden_tests/01_regenie_s1.py -c ${config_file} -i ${input_vcf} --wdir ${wdir} | tee 7_regenie_s1.logs
  """
}