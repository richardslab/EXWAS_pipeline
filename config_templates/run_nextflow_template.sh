#!/bin/bash

######################### TO UPDATE #########################
# specify the nextflow parameters
# output directory for nextflow:
## place for all cache
cache_dir=""

# path to main.nf
nxtflow_script=""
# path to nextflow.config
nxtflow_config=""
##################################################



wdir="${cache_dir}/wdir"
conda_cache_dir="${cache_dir}/conda_cache"
logfile="${cache_dir}/nextflow.log"

mkdir -p $cache_dir
mkdir -p $wdir
mkdir -p $conda_cache_dir

cd $cache_dir
export NXF_LOG_FILE=$logfile
export NXF_WORK=$wdir
export NXF_OFFLINE=true
export NXF_CONDA_CACHEDIR=$conda_cache_dir

# only available in 24.02.0-edge.
# in previous versions, this is in launch dir (which is $cache_dir)
export NXF_CACHE_DIR="${cache_dir}/.nextflow"

# nextflow run $nxtflow_script -c $nxtflow_config -profile standard
nextflow run $nxtflow_script -c $nxtflow_config -profile conda
