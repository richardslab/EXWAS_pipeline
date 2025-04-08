process create_annotation_summaries{
  publishDir("${params.outdir}/ANNOTATE_VARIANTS/04_annotation_summaries",mode:"link")
  errorStrategy "terminate"
  label "big_memory","long"

  input:
    path annotate_vcf_logs
    path annotate_vcf_res
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
  
  output:
    path "5_1_vep_summaries_${ofile_suffix}.sqlite3.db", emit: "annotations_db"
    path "5_1_vep_summaries_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  set -o pipefail
  create_annotation_summaries.py -c ${config_file} -i ${annotation_vcf} | tee "5_1_vep_summaries_${ofile_suffix}.logs"
  """
}