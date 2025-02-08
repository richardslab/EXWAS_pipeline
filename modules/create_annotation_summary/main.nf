process create_annotation_summaries{
  storeDir params.outdir

  input:
    path annotate_vcf
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir
  
  output:
    path "5_1_vep_summaries_${ofile_suffix}.sqlite3.db"
    path "5_1_vep_summaries_${ofile_suffix}.logs",emit: "log"
  
  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/04_1_create_annotation_summaries.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "5_1_vep_summaries_${ofile_suffix}.logs"
  """
}