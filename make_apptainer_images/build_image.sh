#!/bin/bash



module load StdEnv/2023 apptainer/1.2.4


cd /home/richards/kevin.liang2/scratch/exwas_pipeline/results/build_exwas_pipeline_apptainer

cachedir=/home/richards/kevin.liang2/scratch/exwas_pipeline/results/build_exwas_pipeline_apptainer/apptainer_cache
tmpdir=/home/richards/kevin.liang2/scratch/exwas_pipeline/results/build_exwas_pipeline_apptainer/apptainer_tmpdir

def_file=/home/richards/kevin.liang2/scratch/exwas_pipeline/results/build_exwas_pipeline_apptainer/exwas_pipeline_vep105_apptainer.def



export APPTAINER_CACHEDIR=$cachedir
export APPTAINER_TMPDIR=$tmpdir
apptainer build --sandbox --fakeroot test/ $def_file