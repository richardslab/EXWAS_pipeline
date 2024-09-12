/*
  This is the main ExWAS workflow using Regenie.

  Consists of 2 main workflow:
    1. preparing input files for Regenie
    2. Running regenie
  
  Main data analysis are done in python.
  Data inputs are specified in a YAML parameter_file.

  workflow 1: preparing regenie input files
    #' VCF processing
    1. Align input VCF file
    2. Annotate input VCF file with VEP

    #' Regenie inputs
    3. Create mask file for the VCF
    4. Create annotation file for VCF
    5. Create set list file for VCF
*/

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

include {say_hi} from "./modules/test_py"

workflow {
  hi_ch = say_hi('hello')
  hi_ch | view { it }
}



workflow.onComplete {
	log.info (
    workflow.success ? '\nDone!' : '\nPipeline did not complete' )
}