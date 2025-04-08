process create_annotation_file {
  publishDir("${params.outdir}/REGENIE_INPUTS",mode:"link")
  errorStrategy "terminate"
  label "light","instant"

  input:
    val annotation_summaries
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
  
  output:
    path "*/annotations_${ofile_suffix}.txt",emit: "regenie_annotations"
    path "5_2_create_annotation_file_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  set -o pipefail
  create_annotation_files.py -c ${config_file} -i ${annotation_vcf}  --annotation_summaries ${annotation_summaries} | tee "5_2_create_annotation_file_${ofile_suffix}.logs"
  """
}