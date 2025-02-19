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


include {ANNOTATE_VARIANTS} from "../subworkflows/local/Annotation_main.nf"
include {CREATE_REGENIE_INPUT} from "../subworkflows/local/create_regenie_inputs.nf"


workflow EXWAS_INPUT_PREP_PIPELINE {
  take:
    apptainer_img

  main:
    // Generate annotations, which will be input file path, base name tuples
    Channel.fromPath(params.annotation_vcf).map{
        file -> [file,"${file.baseName}"]
    }.set{ annotation_inputs }
    
    // SUBWORKFLOW: Run input validation and variant annotation
    anno_res = ANNOTATE_VARIANTS(
      apptainer_img,
      annotation_inputs
    )

    // // SUBWORKFLOW: create Regenie inputs
    regenie_input_res = CREATE_REGENIE_INPUT(
      anno_res.anno_db,
      anno_res.anno_log,
      annotation_inputs
    )

  emit:
    regenie_input_prep = regenie_input_res.setlist_res.join(
      regenie_input_res.annotate_res
    ).join(
      regenie_input_res.mask_res
    )


}