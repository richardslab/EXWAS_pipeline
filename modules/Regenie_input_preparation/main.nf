process check_yaml_config{
  /** perform sanity checks on the configuration file
  * 
  * will output python script print to stdout and pipe to log file
  */

  storeDir params.outdir
  input:
    val config_file
    val input_vcf
    val wdir

  output:
    path "01_yaml_config_checks.log"

  script:
  """
  python ${baseDir}/modules/Regenie_input_preparation/00_check_config_file.py -c ${config_file} --input_vcf ${input_vcf} --wdir ${wdir} | tee 01_yaml_config_checks.log
  """
}

process align_vcf{
  /** left align VCF
  */
  storeDir params.outdir

  output:
    path "02_vcf_alignment.log"
  
  script:
  """
  """
}