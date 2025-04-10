process align_vcf{
  publishDir("${params.outdir}/ANNOTATE_VARIANTS/02_VCF_alignment",mode:"link")
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
    return baseTime + (task.attempt ?: 1  ) * 1.h
  }

  /** left align VCF to reference
  */
  input:
    path check_yaml_config_logs
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file


  output:
    path "2_vcf_alignment_${ofile_suffix}.log", emit: 'log'
    path "2_bcftool_sitesonly_${ofile_suffix}.vcf.gz", emit: "aligned_vcf_files"
    path "2_bcftool_sitesonly_${ofile_suffix}.vcf.gz.tbi", emit: "aligned_vcf_tabix_index"
  
  script:
  """
  set -o pipefail
  align_format_vcf.py -c ${config_file} -i ${annotation_vcf} | tee "2_vcf_alignment_${ofile_suffix}.log"
  """
}