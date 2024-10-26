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
Prepare and run ExWAS gene burden tests from input VCF files

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
include {check_yaml_config; align_vcf; annotate_vcf; create_mask_files; create_annotation_summaries; create_annotation_file; create_setlist_file} from "./modules/Regenie_input_preparation"

include {run_regenie_s1; run_regenie_s2} from "./modules/Regenie_gene_burden_tests"
workflow testflow{
    process test {
    /*
      Simple process to just echo what it got
    */
    input:
      val x
      tuple val(a),val(b)
      
    output:
      stdout

    script:
    """
    echo ${x} ${a} ${b}
    """
  }
  process test1{
    output:
      stdout
    script:
    """
    echo 'beginning'
    """
  }
  Channel.fromPath(params.step2_exwas_genetic).filter{
    file -> file.name.endsWith(".${params.step2_exwas_genetic_file_type}")
  }.map{
    file -> [file,"${file.baseName}"]
  }.set{ regenie_input }

  test1()
  test("abc",regenie_input,test1.log) | view{it}
}

// This workflow runs ExWAS for 1 VCF
// by using Channel, can run mutliple VCF in parallele (i.e, per chr)
workflow annotation_workflow {
  take:
    input_tuple

  main:
    check_yaml_config(input_tuple,params.config_file,params.outdir)  

    align_vcf(input_tuple,params.config_file,params.outdir,check_yaml_config.out.log)
    
    annotate_vcf(input_tuple,params.config_file,params.outdir,align_vcf.out.log)

    create_mask_files(input_tuple,params.config_file,params.outdir,annotate_vcf.out.log)

    create_annotation_summaries(input_tuple,params.config_file,params.outdir,create_mask_files.out.log)

    create_annotation_file(input_tuple,params.config_file,params.outdir,create_annotation_summaries.out.log)

    create_setlist_file(input_tuple,params.config_file,params.outdir,create_annotation_summaries.out.log)

    emit:
      create_setlist_file.out.log
}

workflow regenie_workflow {
  take:
    annotation_workflow_log
  
  main:
    Channel.fromPath(params.step2_exwas_genetic).filter{
      file -> file.name.endsWith(".${params.step2_exwas_genetic_file_type}")
    }.map{
      file -> [file,"${file.baseName}"]
    }.set{ regenie_input }

    run_regenie_s1(params.config_file,params.outdir)

    run_regenie_s2(regenie_input,params.config_file,params.outdir,params.step2_exwas_genetic,params.step2_exwas_genetic_file_type,params.annotation_vcf,run_regenie_s1.out.log,annotation_workflow_log)
}


workflow {
  // Generate annotations, which will be input file path, base name tuples
  Channel.fromPath(params.annotation_vcf).map{
    file -> [file,"${file.baseName}"]
  }.set{ annotation_inputs }

  annotation_workflow(annotation_inputs)

  // Run Regenie  
  regenie_workflow(annotation_workflow.out)
}


workflow.onComplete {
	log.info (
    workflow.success ? '\nDone!' : '\nPipeline did not complete' )
}
