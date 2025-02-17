process create_setlist_file {
  publishDir("${params.outdir}/REGENIE_INPUTS",mode:"link")
  errorStrategy "terminate"
  label "light","instant"

  input:
    path annotation_summaries
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
  
  output:
    path "6_${ofile_suffix}.setlist", emit: "setlist_files"
    path "6_create_setlist_file_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  set -o pipefail
  dir=\$(dirname ${annotation_summaries})
  create_set_list_file.py -c ${config_file} -i ${annotation_vcf} --annotation_summary_dir \$dir | tee "6_create_setlist_file_${ofile_suffix}.logs"
  """
}