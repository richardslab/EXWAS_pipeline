// This is the workflow to generate summaries of Regenie results

process find_data{
  publishDir("${params.outdir}/SUMMARY",mode:"link")
  errorStrategy "terminate"
  label "light","short"

  input:
    val config
    val regenie_results

  output:
    path "9_found_regenie_results.log", emit: "log"
    path "*/Regenie_Summaries/Phenotypes_results_paths.yaml.gz", emit: "result_paths"

  script:
    """
    set -o pipefail
    find_data.py -c ${config} --idir ${params.outdir}/REGENIE_OUTPUTS/Regenie_S2 | tee 9_found_regenie_results.log
    """
}

process compute_lambda{
  publishDir("${params.outdir}/SUMMARY",mode:"link")
  errorStrategy "terminate"
  label "light","short"

  input:
    val config
    val find_data_pheno_file

  output:
    path "10_lambda_results.log", emit: "log"
    path "*/Regenie_Summaries/Genomic_inflation_factors.tsv.gz", emit: "genomic_inflation_factors"
  
  script:
    """
    set -o pipefail
    compute_lambda.py -c ${config} --res_path ${find_data_pheno_file} | tee 10_lambda_results.log
    """
}

process obtain_assoc_counts{
  publishDir("${params.outdir}/SUMMARY",mode:"link")
  errorStrategy "terminate"
  label "light","short"

  input:
    val config
    val find_data_res

  output:
    path "11_assoc_counts.log",emit: "log"
    path "*/Regenie_Summaries/ExWAS_counts.tsv.gz", emit: "ExWAS_counts"
  
  script:
    """
    set -o pipefail
    association_counts.py -c ${config} --res_path ${find_data_res} | tee 11_assoc_counts.log
    """
}