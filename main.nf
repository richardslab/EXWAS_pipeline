#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    richardslab/ExWAS pipelien
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/richardslab/EXWAS_pipeline
----------------------------------------------------------------------------------------

Runs main workflow
*/

log.info """
    Regenie ExWAS gene burden test pipelien
"""



include {EXWAS_INPUT_PREP_PIPELINE} from "./workflows/exwas_input_prep"
include {BUILD_VEP_IMG} from "./workflows/create_vep_apptainer_img"
include {RUN_REGENIE} from "./workflows/run_regenie"

workflow {
    // WORKFLOW: CREATE VEP APPTAINER IMAGE
    build_res = BUILD_VEP_IMG()
    
    // WORKFLOW: Preparae regenie inputs
    //
    EXWAS_INPUT_PREP_PIPELINE(build_res.apptainer_img)

    // WORKFLOW: Run Regenie
    RUN_REGENIE()

    // WORKFLOW: Summarize Regenie results
    // SUMMARIZE_RESULTS()

}

workflow.onComplete{
  log.info (workflow.success ? '\nCompleted all ExWAS!!' : '\nPipeline did not complete' )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

