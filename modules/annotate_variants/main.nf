process annotate_vcf {
  publishDir("${params.outdir}/ANNOTATE_VARIANTS/03_variant_annotations",mode:"link")
  label "big_memory","short"

  // Retry the job up to 3 times if it fails
  errorStrategy { task.exitStatus in [137, 143, 139, 140] ? 'retry' : 'terminate' }
  maxRetries 3

  // Dynamically increase memory with each retry
  memory {
    def baseMem = 32.GB
    return baseMem * (task.attempt ?: 1)
  }

  // increase runtime if timed out
  time {
    def baseTime = 24.h
    return baseTime + (task.attempt ?: 1  ) * 12.h
  }

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