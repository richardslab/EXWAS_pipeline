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

## Required files:
 * exwas_pipeline.yml: conda environment file to execute the python scripts
   * Specifies required python packages
 * nextflow_template.config: nextflow configuration
   * Specifies:
     * input VCF files
     * output directory
     * ExWAS configurations
 * proj_config_template.yml: ExWAS configuration yaml files:
   * Program executable paths:
     * nextflow >= 23.10.0
     * bcftools
     * tabix
     * plink
     * plink2
     * vep apptainer image
     * vep cache directory
     * regenie executable
     * python script and heper files
   * ExWAS parameters
     * Mask definitions
