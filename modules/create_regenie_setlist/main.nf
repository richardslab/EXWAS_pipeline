process create_setlist_file {
  publishDir("${params.outdir}/REGENIE_INPUTS",mode:"link")

  // Retry the job up to 3 times if it fails
  errorStrategy { task.exitStatus in [137, 143, 139, 140] ? 'retry' : 'terminate' }
  maxRetries 3

  // Dynamically increase memory with each retry
  memory {
    def baseMem = 15.GB
    return baseMem * (task.attempt ?: 1)
  }

  // increase runtime if timed out
  time {
    def baseTime = 1.h
    return baseTime + (task.attempt ?: 1  ) * 2.h
  }

  input:
    val annotation_summaries
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
  
  output:
    path "6_${ofile_suffix}.setlist", emit: "setlist_files"
    path "6_create_setlist_file_${ofile_suffix}.logs", emit: "log"
  
  script:
  """
  set -o pipefail
  create_set_list_file.py -c ${config_file} -i ${annotation_vcf} --idir ${params.outdir}/ANNOTATE_VARIANTS/04_annotation_summaries | tee "6_create_setlist_file_${ofile_suffix}.logs"
  """
}