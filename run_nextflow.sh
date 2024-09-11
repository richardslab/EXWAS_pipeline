#!/bin/zsh
source ~/.zshrc
module load StdEnv/2023 nextflow/23.10.0

######################### TO UPDATE #########################
# specify the nextflow cache and log directories
# these are environment variables of Nextflow so have to specify in launching environment (i.e., here)
cache_dir="/scratch/richards/kevin.liang2/exwas_pipeline/results/nextflow_out"
nxtflow_script="/home/richards/kevin.liang2/scratch/exwas_pipeline/src/main.nf"
nxtflow_config="/home/richards/kevin.liang2/scratch/exwas_pipeline/src/nextflow.config"

##################################################



wdir="${cache_dir}/wdir"
logfile="${cache_dir}/nextflow.log"

mkdir -p $cache_dir
mkdir -p $wdir
cd $cache_dir
export NXF_LOG_FILE=$logfile
export NXF_WORK=$wdir
export NXF_OFFLINE=true
# only available in 24.02.0-edge
export NXF_CACHE_DIR="${cache_dir}/.nextflow"

nextflow run $nxtflow_script -c $nxtflow_config -profile conda
