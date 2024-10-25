// This is the workflow to generate summaries of Regenie results


process find_data{
  input:
    val config
    val wdir

  output:
    path "9_found_regenie_results.log", emit: log

  script:
    """
    set -o pipefail

    python -u ${baseDir}/modules/Regenie_results_summaries/00_find_data.py -c ${config} --wdir ${wdir} | tee 9_found_regenie_results.log
    """
}

process compute_lambda{
  input:
    val config
    val wdir
    val find_data_logs

  output:
    path "10_lambda_results.log"
  
  script:
    """
    set -o pipefail

    python -u ${baseDir}/modules/Regenie_results_summaries/01_compute_lambda.py -c ${config} --wdir ${wdir} | tee 10_lambda_results.log
    """
}