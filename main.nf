#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    richardslab/ExWAS pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/richardslab/EXWAS_pipeline
----------------------------------------------------------------------------------------

Runs main workflow
*/

log.info """
    Regenie ExWAS gene burden test pipeline
"""



include {EXWAS_INPUT_PREP_PIPELINE} from "./workflows/exwas_input_prep"
include {BUILD_VEP_IMG} from "./workflows/create_vep_apptainer_img"
include {RUN_REGENIE} from "./workflows/run_regenie"

include {SUMMARIZE_RESULTS} from "./workflows/summarize_results"

workflow {
    // WORKFLOW: CREATE VEP APPTAINER IMAGE
    build_res = BUILD_VEP_IMG()
    
    // WORKFLOW: Preparae regenie inputs
    exwas_input_prep_res = EXWAS_INPUT_PREP_PIPELINE(
        build_res.apptainer_img
    )

    // WORKFLOW: Run Regenie
    regenie_res = RUN_REGENIE(
        exwas_input_prep_res.regenie_input_prep_mask.collect().map{it.sort()},
        exwas_input_prep_res.regenie_input_prep_setlist.collect().map{it.sort()},
        exwas_input_prep_res.regenie_input_prep_annotation.collect().map{it.sort()}
    )

    // WORKFLOW: Summarize Regenie results
    SUMMARIZE_RESULTS(
        regenie_res.regenie_s2.collect().map{it.sort()}
    )

}

workflow.onComplete{
  log.info (workflow.success ? '\nCompleted all ExWAS!!' : '\nPipeline did not complete' )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

