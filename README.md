# EXWAS_pipeline

## running it (for now)
1. put the conda environment file (exwas_pipeline.yml) in this directory at the same level as README.md
     * will be provided probably once it is done
3. Fill in nextflow_template.config and proj_config_template.yml
4. run:
```
nextflow run <path>/main.nf -c <path>/nextflow_template.config -profile conda
```
OR edit run_nextflow_template.sh with proper in/out directories for nextflow. then run
```
<path>/run_nextflow_template.sh
```

Making conda environment on first run will take some time. As long as the conda cache dir is not deleted, the environment will not be made again.

## Required files (examples in config_templates):
 * exwas_pipeline.yml: conda environment file to execute the python scripts
 * nextflow_template.config: nextflow configuration
 * proj_config_template.yml: ExWAS configuration yaml files:
## program requirements:
 * nextflow >= 23.10.0
 * python 3.10.9
 * apptainer >= 1.2.4
 * bcftools
 * tabix
 * plink
 * plink2
 * vep apptainer image (probably will be provided, or make from definition file)
 * regenie
