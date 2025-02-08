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

workflow {
    // WORKFLOW: Run main workflow
    //
    EXWAS_PIPELINE()
}

workflow.onComplete{
  log.info (workflow.success ? '\nCompleted all ExWAS!!' : '\nPipeline did not complete' )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

