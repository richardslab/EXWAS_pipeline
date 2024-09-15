// This is the base file name produced by the scripts
// required so nextflow can check all outputs are produced
def vcf_file = new File(params.input_vcf)
def vcf_file_name = vcf_file.getName()

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
  # this is required because even if python script crashes, the tee will work so the final error code is fine and nextflow doesn't know.

  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/00_check_config_file.py -c ${config_file} --input_vcf ${input_vcf} --wdir ${wdir} | tee 01_yaml_config_checks.log
  """
}



process align_vcf{
  /** left align VCF to reference
  */
  storeDir params.outdir
  input:
    val config_file
    val input_vcf
    val wdir

  output:
    
    path "02_vcf_alignment.log"
    path "1_${vcf_file_name}_bcftool_variant_only.set_ids.no_genotypes.vcf.gz"
    path "1_${vcf_file_name}_bcftool_variant_only.set_ids.no_genotypes.vcf.gz.tbi"
  
  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/01_align_format_vcf.py -c ${config_file} -i ${input_vcf} --wdir ${wdir} | tee 02_vcf_alignment.log
  """
}

process annotate_vcf {
  storeDir params.outdir
  input:
    val config_file
    val input_vcf
    val wdir

  output:
    path "03_vep_annotation.logs"
    // path "2_${vcf_file_name}_vcf_final_annotation.txt"

  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/02_annotate_with_vep.py -c ${config_file} -i ${input_vcf} --wdir ${wdir} | tee 03_vep_annotation.logs
  """

}