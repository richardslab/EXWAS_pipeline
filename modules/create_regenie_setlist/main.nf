process create_setlist_file {
  storeDir params.outdir

  input:
    path create_annotation_summaries_logs
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir
  
  output:
    path "6_${ofile_suffix}.setlist"
    path "6_create_setlist_file_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/05_create_set_list_file.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "6_create_setlist_file_${ofile_suffix}.logs"
  """
}