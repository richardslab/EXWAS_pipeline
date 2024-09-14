# EXWAS_pipeline

## Required files:
 * conda environment file
   * Specifies required python packages
 * nextflow.config file
   * Specifies:
     * input VCF files
     * output directory
     * ExWAS configurations
 * ExWAS configuration yaml files:
   * Program executable paths:
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