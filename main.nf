#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    richardslab/ExWAS pipelien
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/richardslab/EXWAS_pipeline
----------------------------------------------------------------------------------------

Runs main workflow
*/



include {EXWAS_PIPELINE} from "./workflows/exwas"
include {BUILD_VEP_IMG} from "./workflows/create_vep_apptainer_img.nf"

workflow {
    // WORKFLOW: CREATE VEP APPTAINER IMAGE
    build_res = BUILD_VEP_IMG()
    
    // WORKFLOW: Run main workflow
    //
    EXWAS_PIPELINE(build_res.apptainer_img)

}

workflow.onComplete{
  log.info (workflow.success ? '\nCompleted all ExWAS!!' : '\nPipeline did not complete' )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

