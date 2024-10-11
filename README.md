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

## Usage notes
  * Specified within nextflow_template.config:
    * VCF for generating annotation files are specified separately from the input to run Regenie. Sites-only VCF files can be used to generate annotation, as it is smaller file size and faster to run. As long as the variants are the same, should not need to regenerate the annotations for each ExWAS (e.g., same set of annotations for all males and females, stratified analyses, etc)
      * if these VCF files have genetic data, a sites only vcf file will be created in any case using bcftools
    * Annotation files generated from VCF files will be matched by wildcard character if specified or assumes a 1-1 matching
        * e.g., if annotation file has name: **Sites_only_VCF_chr\*.vcf** and the regenie input file has name **Another_file_chr\*.pgen**, then will match based on whatever is specified by the character in the '*' position
        * e.g., if annotation file has name: **Sites_only_allchr.vcf** and the regenie input is **Another_file_allchr.bgen**, then it will generate only 1 annotation file and assumes it matches to **Another_file_allchr.bgen**
        * the wildcard character can stand-in for 1 or more alphanumeric symbols.
    * Regenie input ExWAS file is passed directly to Regenie, so specify as it should based on [Regenie documentation](https://rgcgithub.github.io/regenie/options/) (e.g., pgen/bfile have no extensions)

## Configuration files
  * exwas_pipeline.yml: conda environment file to execute the python scripts
  * nextflow_template.config: nextflow configuration
  * proj_config_template.yml: ExWAS configuration yaml files:
## program requirements (paths to be specified in proj_config_template.yml):
  * nextflow >= 23.10.0
  * python 3.10.9
  * Other required programs (plink, tabix, etc) are listed in proj_config_template.yml
