process annotate_vcf {
  publishDir("${params.outdir}/ANNOTATE_VARIANTS/03_variant_annotations",mode:"link")
  errorStrategy "terminate"
  label "big_memory","long"

  input:
    path align_vcf_logs
    path align_vcfs
    path align_vcf_tabix
    val  apptainer_img
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file

  output:
    path "3_vep_annotation_${ofile_suffix}.logs", emit: 'log'
    path "3_annotation_results_${ofile_suffix}.txt.gz", emit: "var_annotations"

  script:
  """
  set -o pipefail
  # Because it is calling apptainer img with mounts, pass in 
  # directory of published file to access directly.
  # things in workdir were symlinks and will not work.
  annotate_with_vep.py -c ${config_file} -i ${annotation_vcf} --vep_img ${apptainer_img} --input_dir ${params.outdir}/ANNOTATE_VARIANTS/02_VCF_alignment | tee "3_vep_annotation_${ofile_suffix}.logs"
  """
}