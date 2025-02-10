/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Variant annotation VEP
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This subworkflow annootate variants using VEP with user defined plugins.
*/

// Uses DSL 2
nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check runtime environment.
Conda should not be included
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Check input files.
If wildcards are used for VCF file for VEP annotations, wildcards should be present in the 
genetic files for burden testing also
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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
              Variant annotation Subworkflow
==================================================
Annotate variant using VEP
==================================================
"""


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LOAD IN REQUIRED MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include {check_yaml_config} from "../../modules/validate_config"
include {align_vcf} from "../../modules/align_vcf"
include {annotate_vcf} from "../../modules/annotate_variants"
include {create_annotation_summaries} from "../../modules/create_annotation_summary"
include{build_vep_apptainer_img} from "../../modules/create_vep_apptainer_img"

workflow ANNOTATE_VARIANTS {
  take:
    apptainer_img
    each_input

  main:
    check_res = check_yaml_config(
      each_input,
      params.config_file
    )

    // align_vcf_res = align_vcf(
    //   check_res.log.collect(),
    //   each_input,
    //   params.config_file
    // )
    
    // annotate_res = annotate_vcf(
    //   align_vcf_res.log.collect(),
    //   align_vcf_res.aligned_vcf_files.collect(),
    //   align_vcf_res.aligned_vcf_tabix_index.collect(),
    //   apptainer_img,
    //   each_input,
    //   params.config_file
    // )

    // annotate_summary_res = create_annotation_summaries(
    //   annotate_res.log.collect(),
    //   annotate_res.var_annotations.collect(),
    //   each_input,
    //   params.config_file
    // )
}