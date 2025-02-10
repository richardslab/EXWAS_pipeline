process annotate_vcf {
  publishDir("${params.outdir}/ANNOTATE_VARIANTS/03_variant_annotations",mode:"link")
  errorStrategy "terminate"
  label "big_memory","long"

  input:
    path align_vcf_logs
    path align_vcfs
    path align_vcf_tabix
    val apptainer_img
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file

  output:
    path "3_vep_annotation_${ofile_suffix}.logs", emit: 'log'
    path "3_annotation_results_${ofile_suffix}.txt.gz", emit: "var_annotations"

  script:
  """
  set -o pipefail
  annotate_with_vep.py -c ${config_file} -i ${annotation_vcf} --vep_img ${apptainer_img} | tee "3_vep_annotation_${ofile_suffix}.logs"
  """
}