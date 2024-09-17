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



def checkRuntimeEnvironment(){
  /**
  * Check the runtime environment for nextflow before any process
  * 1. It should not have any conda environment activated because it messes up the paths
  */
  python_vrs = "which python".execute().text
  conda_prefix = System.getenv("CONDA_PREFIX")
  return [python_vrs,conda_prefix]
}
def runtimeEnv = checkRuntimeEnvironment()
assert runtimeEnv[1] == null : "Deactivate conda environment first. Founrd ${runtimeEnv[1]}"





log.info """
              Regenie ExWAS pipeline
==================================================
Prepare and run ExWAS gene burden tests from input VCF files

  # input data
  Input vcf file: ${params.input_vcf}

  # output data
  output directory: ${params.outdir}

  # pipeline configurations
  conda environment: ${baseDir}/exwas_pipeline_conda_env.yml
  pipeline configuration: ${params.config_file}

"""

// Regenie input processing
include {check_yaml_config; align_vcf; annotate_vcf; create_mask_files; create_annotation_summaries; create_annotation_file; create_setlist_file} from "./modules/Regenie_input_preparation"

include {run_regenie_s1} from "./modules/Regenie_gene_burden_tests"

workflow regenie_workflow {

  check_yaml_config(params.config_file,params.input_vcf,params.outdir)  

  align_vcf(params.config_file,params.input_vcf,params.outdir,check_yaml_config.out.log)
  
  annotate_vcf(params.config_file,params.input_vcf,params.outdir,align_vcf.out.log)

  create_mask_files(params.config_file,params.input_vcf,params.outdir,annotate_vcf.out.log)

  create_annotation_summaries(params.config_file,params.input_vcf,params.outdir,create_mask_files.out.log)

  create_annotation_file(params.config_file,params.input_vcf,params.outdir,create_annotation_summaries.out.log)

  create_setlist_file(params.config_file,params.input_vcf,params.outdir,create_annotation_summaries.out.log)
  
  // run_regenie_s1(params.config_file,params.input_vcf,params.outdir,create_setlist_file.out.log)
}

workflow {
  regenie_workflow()
}



workflow.onComplete {
	log.info (
    workflow.success ? '\nDone!' : '\nPipeline did not complete' )
}