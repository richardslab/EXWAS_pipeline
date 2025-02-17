// This is the workflow to generate summaries of Regenie results

process find_data{
  publishDir("${params.outdir}/SUMMARY",mode:"link")
  errorStrategy "terminate"
  label "light","short"

  input:
    val config
    val wdir

  output:
    path "9_found_regenie_results.log", emit: log
    path "Phenotypes_results_paths.yaml.gz", emit: result_paths

  script:
    """
    set -o pipefail
    find_data.py -c ${config} --wdir ${wdir} | tee 9_found_regenie_results.log
    """
}

process compute_lambda{
  publishDir("${params.outdir}/SUMMARY",mode:"link")
  errorStrategy "terminate"
  label "light","short"

  input:
    val config
    val wdir
    val find_data_logs

  output:
    path "10_lambda_results.log", emit: log
  
  script:
    """
    set -o pipefail
    compute_lambda.py -c ${config} --wdir ${wdir} | tee 10_lambda_results.log
    """
}

process obtain_assoc_counts{
  publishDir("${params.outdir}/SUMMARY",mode:"link")
  errorStrategy "terminate"
  label "light","short"

  input:
    val config
    val wdir
    val find_data_logs

  output:
    path "11_assoc_counts.log"
  
  script:
    """
    set -o pipefail
    association_counts.py -c ${config} --wdir ${wdir} | tee 11_assoc_counts.log
    """
}