#!/bin/zsh
source ~/.zshrc


# specify the nextflow cache and log directories
# these are environment variables of Nextflow so have to specify in launching environment (i.e., here)
export NXF_LOG_FILE="/Users/kevinliang/Desktop/work/working/exwas_pipelines/nextflow_cache"
export NXF_CACHE_DIR="/Users/kevinliang/Desktop/work/working/exwas_pipelines/nextflow_cache/.nextflow"
export NXF_WORK="/Users/kevinliang/Desktop/work/working/exwas_pipelines/nextflow_cache"

export NXF_OFFLINE=true

nextflow run /Users/kevinliang/Desktop/work/working/exwas_pipelines/src/test_nested_workflow/main.nf