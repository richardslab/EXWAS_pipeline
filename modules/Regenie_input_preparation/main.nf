// This is the base file name produced by the scripts
// required so nextflow can check all outputs are produced

process check_yaml_config{
  /** perform sanity checks on the configuration file
  * 
  * will output python script print to stdout and pipe to log file
  */
  storeDir params.outdir
  

  input:
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir
    

  output:
    path "1_yaml_config_checks_${ofile_suffix}.log", emit: 'log'

  script:
  """
  # this is required because even if python script crashes, the tee will work so the final error code is fine and nextflow doesn't know.

  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/00_check_config_file.py -c ${config_file} --input_vcf ${annotation_vcf} --wdir ${wdir} | tee "1_yaml_config_checks_${ofile_suffix}.log"
  """
}



process align_vcf{
  /** left align VCF to reference
  */
  storeDir params.outdir
  input:
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir
    path "1_yaml_config_checks_${ofile_suffix}.log"


  output:
    path "2_vcf_alignment_${ofile_suffix}.log", emit: 'log'
    path "2_bcftool_sitesonly_${ofile_suffix}.vcf.gz"
    path "2_bcftool_sitesonly_${ofile_suffix}.vcf.gz.tbi"
  
  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/01_align_format_vcf.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "2_vcf_alignment_${ofile_suffix}.log"
  """
}

process annotate_vcf {
  storeDir params.outdir
  input:
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir
    path "2_vcf_alignment_${ofile_suffix}.log"

  output:
    path "3_vep_annotation_${ofile_suffix}.logs", emit: 'log'
    path "3_annotation_results_${ofile_suffix}.txt"

  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/02_annotate_with_vep.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "3_vep_annotation_${ofile_suffix}.logs"
  """
}

process create_mask_files {
  storeDir params.outdir
  input: 
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir
    path "3_vep_annotation_${ofile_suffix}.log"
  
  output:
    // wildcard for study names
    // figure out if it can specify numbers
    path "*/masks_${ofile_suffix}.txt" 
    path "4_create_masks_${ofile_suffix}.logs", emit: "log"
  
  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/03_create_mask_files.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "4_create_masks_${ofile_suffix}.logs"
  """
}

process create_annotation_summaries{
  storeDir params.outdir

  input:
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir
    path "4_create_masks_${ofile_suffix}.log"
  
  output:
    path "5_1_vep_summaries_${ofile_suffix}.sqlite3.db"
    path "5_1_vep_summaries_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/04_1_create_annotation_summaries.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "5_1_create_masks_${ofile_suffix}.logs"
  """
}

process create_annotation_file {
  storeDir params.outdir

  input:
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir
    path "5_1_vep_summaries_${ofile_suffix}.log"
  
  output:
    path "*/annotations_${ofile_suffix}.txt"
    path "5_2_create_annotation_file_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/04_2_create_annotation_files.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "5_2_create_annotation_file_${ofile_suffix}.logs"
  """

}

process create_setlist_file {
  storeDir params.outdir

  input:
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir
    path "5_1_vep_summaries_${ofile_suffix}.log"
  
  output:
    path "6_${ofile_suffix}.setlist"
    path "6_create_setlist_file_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/05_create_set_list_file.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "6_create_setlist_file_${ofile_suffix}.logs"
  """
}