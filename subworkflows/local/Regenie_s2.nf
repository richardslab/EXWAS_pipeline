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

log.info """
      Regenie gene burden test subworkflow
==================================================
- Run Regenie Step 2 per study defined by user
==================================================
"""

// Regenie input processing
include {run_regenie_s2} from "../../modules/run_regenie_s2"
workflow RUN_REGENIE_S2 {
  take:
    regenie_s1
    regenie_inputs_mask
    regenie_inputs_setlist
    regenie_inputs_annotation

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

    
    regenie_s2_res = run_regenie_s2(
      regenie_input,
      regenie_inputs_mask,
      regenie_inputs_setlist,
      regenie_inputs_annotation,
      params.config_file,
      regenie_s1,
      params.step2_exwas_genetic,
      params.step2_exwas_genetic_file_type,
      params.annotation_vcf
    )
    
  emit:
    regenie_s2 = regenie_s2_res.regenie_S2_results
    regenie_s2_logs = regenie_s2_res.regenie_logs

}
