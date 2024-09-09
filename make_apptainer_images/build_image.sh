#!/bin/bash

cd ~

def_file=/tmp/lima/exwas_pipeline_vep105_apptainer.def



apptainer build --fakeroot vep_apptainer $def_file