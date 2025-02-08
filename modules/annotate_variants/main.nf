process annotate_vcf {
  storeDir params.outdir
  input:
    path align_vcf_logs
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir

  output:
    path "3_vep_annotation_${ofile_suffix}.logs", emit: 'log'
    path "3_annotation_results_${ofile_suffix}.txt.gz"

  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/02_annotate_with_vep.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "3_vep_annotation_${ofile_suffix}.logs"
  """
}