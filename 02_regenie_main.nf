/**
*  This is the main ExWAS workflow using Regenie.
*
*  Consists of 2 parts:
*    1. preparing input files for Regenie
*    2. Running regenie
*  
*  Main data analysis are done in python.
*  Data inputs are specified in a YAML parameter_file.
*
*  1: preparing regenie input files
*    #' VCF processing
*    1. Align input VCF file
*    2. Annotate input VCF file with VEP
*
*    #' Regenie inputs
*    3. Create mask file for the VCF
*    4. Create annotation file for VCF
*    5. Create set list file for VCF
*  2. Run regenie burden testing
*/

// Uses DSL 2
nextflow.enable.dsl = 2

def checkRuntimeEnvironment(){
  /*
   Check the runtime environment for nextflow before any process
   1. It should not have any conda environment activated because it messes up the paths
  */
  python_vrs = "which python".execute().text
  conda_prefix = System.getenv("CONDA_PREFIX")
  return [python_vrs,conda_prefix]
}
def runtimeEnv = checkRuntimeEnvironment()
assert runtimeEnv[1] == null : "Deactivate conda environment first. Found ${runtimeEnv[1]}"

def check_input_exwas_files(){
  /*
    Determine if the wild card character is present in the input files
  */
  pattern = /\*/
  vcf_wildcard = params.annotation_vcf =~ pattern
  exwas_wildcard = params.step2_exwas_genetic =~ pattern
  return [vcf_wildcard,exwas_wildcard]
}
def wildcard_found = check_input_exwas_files()

assert (wildcard_found[0] && wildcard_found[1]) || (! wildcard_found[0] && ! wildcard_found[1]) : "Wildcard notation not consistent between input VCF and ExWAS VCF"
if (wildcard_found[0] && wildcard_found[1]){
  log_str = "Will match annotation VCF with ExWAS input file based on wild card location"
}else{
  log_str = "Assumes 1-1 matchign between VCF and ExWAS input"
}

log.info """
              Regenie ExWAS pipeline
==================================================
Run ExWAS burden testing with Regenie

  # input data
  Input vcf file for generating annotations: ${params.annotation_vcf}
  File for ExWAS: ${params.step2_exwas_genetic}
    ExWAS file type: ${params.step2_exwas_genetic_file_type}
    ${log_str}
  # output data
  output directory: ${params.outdir}

  # pipeline configurations
  conda environment: ${baseDir}/exwas_pipeline_conda_env.yml
  pipeline configuration: ${params.config_file}
==================================================
"""

// Regenie input processing
include {run_regenie_s1; run_regenie_s2} from "./modules/Regenie_gene_burden_tests"

workflow regenie_workflow {
  main:
    if (wildcard_found[1]){
      Channel.fromPath(params.step2_exwas_genetic).filter{
        file -> file.name.endsWith(".${params.step2_exwas_genetic_file_type}")
      }.map{
        file -> [file,"${file.baseName}"]
      }.set{ regenie_input }
    }else{
      filename = file(params.step2_exwas_genetic + "." + params.step2_exwas_genetic_file_type)
      Channel.of(filename).map{
        file -> [file,"${file.baseName}"]
      }.set{ regenie_input }
    }

    run_regenie_s1(params.config_file,params.outdir)
    
    run_regenie_s2(
      regenie_input,
      params.config_file,
      params.outdir,
      params.step2_exwas_genetic,
      params.step2_exwas_genetic_file_type,
      params.annotation_vcf,
      run_regenie_s1.out.log
    )
}

workflow {
  regenie_workflow()
}


workflow.onComplete{
  log.info (workflow.success ? '\nDone Regenie!' : '\nPipeline did not complete' )
}


