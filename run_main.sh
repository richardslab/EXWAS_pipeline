#!/bin/zsh
source ~/.zshrc


# specify the nextflow cache and log directories
export NXF_LOG_FILE="/Users/kevinliang/Desktop/work/working/exwas_pipelines/nextflow_cache"
export NXF_CACHE_DIR="/Users/kevinliang/Desktop/work/working/exwas_pipelines/nextflow_cache/.nextflow"
export NXF_WORK="/Users/kevinliang/Desktop/work/working/exwas_pipelines/nextflow_cache"

export NXF_OFFLINE=true

nextflow run /Users/kevinliang/Desktop/work/working/exwas_pipelines/src/main.nf -c /Users/kevinliang/Desktop/work/working/exwas_pipelines/src/nextflow.config