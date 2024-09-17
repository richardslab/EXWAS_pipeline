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
    path "1_yaml_config_checks.log", emit: 'log'

  script:
  """
  # this is required because even if python script crashes, the tee will work so the final error code is fine and nextflow doesn't know.

  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/00_check_config_file.py -c ${config_file} --input_vcf ${input_vcf} --wdir ${wdir} | tee 1_yaml_config_checks.log
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
    path "1_yaml_config_checks.log"

  output:
    
    path "2_vcf_alignment.log", emit: 'log'
    path "2_${vcf_file_name}_bcftool_variant_only.set_ids.no_genotypes.vcf.gz"
    path "2_${vcf_file_name}_bcftool_variant_only.set_ids.no_genotypes.vcf.gz.tbi"
  
  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/01_align_format_vcf.py -c ${config_file} -i ${input_vcf} --wdir ${wdir} | tee 2_vcf_alignment.log
  """
}

process annotate_vcf {
  storeDir params.outdir
  input:
    val config_file
    val input_vcf
    val wdir
    path "2_vcf_alignment.log"

  output:
    path "3_vep_annotation.logs", emit: 'log'
    path "3_${vcf_file_name}_vcf_final_annotation.txt"

  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/02_annotate_with_vep.py -c ${config_file} -i ${input_vcf} --wdir ${wdir} | tee 3_vep_annotation.logs
  """
}

process create_mask_files {
  storeDir params.outdir
  input: 
    val config_file
    val input_vcf
    val wdir
    path "3_vep_annotation.logs"
  
  output:
    // wildcard for study names
    // figure out if it can specify numbers
    path "*/${vcf_file_name}_masks.txt" 
    path "4_create_masks.logs", emit: "log"
  
  script:
  """
  python -u ${baseDir}/modules/Regenie_input_preparation/03_create_mask_files.py -c ${config_file} -i ${input_vcf} --wdir ${wdir} | tee 4_create_masks.logs
  """
}

process create_annotation_summaries{
  storeDir params.outdir

  input:
    val config_file
    val input_vcf
    val wdir
    path "4_create_masks.logs"
  
  output:
    path "5_1_${vcf_file_name}_vep_summaries.sqlite3.db"
    path "5_1_create_masks.logs",emit: "log"
  
  script:
  """
  python -u ${baseDir}/modules/Regenie_input_preparation/04_1_create_annotation_summaries.py -c ${config_file} -i ${input_vcf} --wdir ${wdir} | tee 5_1_create_masks.logs
  """
}

process create_annotation_file {
  storeDir params.outdir

  input:
    val config_file
    val input_vcf
    val wdir
    path "5_1_create_masks.logs"
  
  output:
    path "*/*_annotations.txt"
    path "5_2_create_annotation_file.logs",emit: "log"
  
  script:
  """
  python -u ${baseDir}/modules/Regenie_input_preparation/04_2_create_annotation_files.py -c ${config_file} -i ${input_vcf} --wdir ${wdir} | tee 5_2_create_annotation_file.logs
  """

}

process create_setlist_file {
  storeDir params.outdir

  input:
    val config_file
    val input_vcf
    val wdir
    path "5_2_create_annotation_file.logs"
  
  output:
    path "6_${vcf_file_name}.setlist"
    path "6_create_setlist_file.logs",emit: "log"
  
  script:
  """
  python -u ${baseDir}/modules/Regenie_input_preparation/05_create_set_list_file.py -c ${config_file} -i ${input_vcf} --wdir ${wdir} | tee 6_create_setlist_file.logs
  """

}