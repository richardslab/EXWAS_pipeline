process create_annotation_file {
  storeDir params.outdir

  input:
    path create_annotation_summaries_logs
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir
  
  output:
    path "*/annotations_${ofile_suffix}.txt"
    path "5_2_create_annotation_file_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/04_2_create_annotation_files.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "5_2_create_annotation_file_${ofile_suffix}.logs"
  """
}