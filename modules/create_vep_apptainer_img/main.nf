process build_vep_apptainer_img{
  publishDir("${params.outdir}/00_VEP_APPTAINER",mode:"link")
  errorStrategy "terminate"
  label "light","short"
  /** left align VCF to reference
  */
  input:
    val config_file
    val def_file

  output:
    path "vep_apptainer.sif", emit: 'apptainer_img'
    path "00_build_vep_apptainer_img.log", emit: 'log'
  
  script:
  """
  set -o pipefail
  create_vep_apptainer_img.py -c ${config_file} --vep_apptainer_def_file ${def_file} | tee "00_build_vep_apptainer_img.log"
  """
}