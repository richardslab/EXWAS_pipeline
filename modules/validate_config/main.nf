process check_yaml_config{
  publishDir("${params.outdir}/ANNOTATE_VARIANTS/01_config_validations",mode:"link")
  errorStrategy "terminate"
  label "light","instant"
  /** perform sanity checks on the configuration file
  * 
  * will output python script print to stdout and pipe to log file
  */
  input:
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    

  output:
    path "1_yaml_config_checks_${ofile_suffix}.log", emit: 'log'

  script:
  """
  # this is required because even if python script crashes, the tee will work so the final error code is fine and nextflow doesn't know.

  set -o pipefail
  check_config_file.py -c ${config_file} --input_vcf ${annotation_vcf} | tee "1_yaml_config_checks_${ofile_suffix}.log"
  """
}