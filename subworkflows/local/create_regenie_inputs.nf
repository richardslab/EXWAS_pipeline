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
              Regenie ExWAS pipeline
==================================================
Annotate and prepare data for ExWAS burden testing with Regenie

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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LOAD IN REQUIRED MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include {create_annotation_summaries} from "./modules/create_annotation_summary"
include {create_annotation_file} from "./modules/create_regenie_annotations"
include {create_setlist_file} from "./modules/create_regenie_setlist"

workflow annotation_workflow {
  take:
    each_input
  
  main:
    // mask_res = create_mask_files(annotate_res.log.collect(),each_input,params.config_file,params.outdir)

    // annotate_file_res = create_annotation_file(annotate_summary_res.log.collect(),each_input,params.config_file,params.outdir)

    // setlist_res = create_setlist_file(annotate_file_res.log.collect(),each_input,params.config_file,params.outdir)

  publish:
      check_res >> "yaml_validations"
      // align_vcf >> "vcf_alignment_results"
      // annotate_res >> "VEP_annotation_results"
      // mask_res >> "Regenie_mask_input_files"
      // annotate_summary_res >> "VEP_annotation_summaries"
      // annotate_file_res >> "Regenie_annotation_inputs_logs"
      // setlist_res >> "Regenie_setlist_input_logs"
}

workflow {
  // Generate annotations, which will be input file path, base name tuples
  Channel.fromPath(params.annotation_vcf).map{
    file -> [file,"${file.baseName}"]
  }.set{ annotation_inputs }
  annotation_workflow(annotation_inputs)
}
workflow.onComplete{
  log.info (workflow.success ? '\nVariant annotation and summary completed!' : '\nPipeline did not complete' )
}


