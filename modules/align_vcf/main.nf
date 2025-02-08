process align_vcf{
  publishDir("yaml_validations",mode:"link")
  /** left align VCF to reference
  */
  storeDir params.outdir
  input:
    path check_yaml_config_logs
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir


  output:
    path "2_vcf_alignment_${ofile_suffix}.log", emit: 'log'
    path "2_bcftool_sitesonly_${ofile_suffix}.vcf.gz", emit: "2_aligned_vcf_files"
    path "2_bcftool_sitesonly_${ofile_suffix}.vcf.gz.tbi", emit: "2_aligned_vcf_tabix_index"
  
  script:
  """
  set -o pipefail
  align_format_vcf.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "2_vcf_alignment_${ofile_suffix}.log"
  """
}