/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
  log_str = "Assumes 1-1 matching between VCF and ExWAS input"
}

include {RUN_REGENIE_S1} from "../subworkflows/local/Regenie_s1"
include {RUN_REGENIE_S2} from "../subworkflows/local/Regenie_s2"


workflow RUN_REGENIE {

  // SUBWORKFLOW: Run step 1 of regenie
  regenie_s1_res = RUN_REGENIE_S1()

  // SUBWORKFLOW: Run step 2 of regenie
  RUN_REGENIE_S2(regenie_s1_res.s1_output)
  
}