# EXWAS_pipeline

## running it (for now)
1. put the conda environment file (exwas_pipeline.yml) at the same level as this main.nf (i.e., where this README.md is)
2. Fill in nextflow_template.config and proj_config_template.yml
3. run:
```
nextflow run <path>/main.nf -c <path>/nextflow_template.config -profile conda
```
  * Nextflow creates caches for runs and stuff, so doing this will store it in the default location.
    
OR edit run_nextflow_template.sh with proper in/out directories for nextflow. then run
```
<path>/run_nextflow_template.sh
```
  * There are paths where you can specify where Nextflow cache should be so you can keep track.
    

## Pipeline notes:
### General
  * Making conda environment on first run will take some time. As long as the conda cache is not deleted, the environment will not be made again.
  * Step 1 expects 1 set of plink files for genotyped data
  
### Annotation
  * VCF for generating annotation files are specified separately from the input to run Regenie. 
      * This way, sites-only VCF files can be used to generate annotation, as it is smaller file size and faster to run.
      * As long as the variants are the same, should not need to regenerate the annotations for each ExWAS (e.g., same set of annotations for all males and females, stratified analyses, etc) (I think...)
    * If VCF contains genetic data, a sites-only VCF file will be created by the --drop-genotypes flag in bcftools. If the VCF is a sites-only VCF file, then this flag will simply not have any effect (I think...)
  * Using VEP
    * The VEP image have no plugins and none of the cache files required to run any plugins. It only has vep installed.
      * Have to download everything, then specify these location in *proj_config_template.yml*
  * The plugins that are parsed right now:
     * IMPACT (HC vs LC)
     * LoFtee
     * CADD
     * dbNSFP
       * alphamissense_pred
       * EVE_Class25_pred
       * LRT_pred
       * MutationTaster_pred
       * Polyphen2_HDIV_pred
       * Polyphen2_HVAR_pred
       * SIFT4G_pred
       * SIFT_pred
      
### ExWAS with Regenie
  * Step 1 of Regenie is done only once and will be used for all 'study' specified in the *proj_config_template.yml*
  * Step 2 of Regenie will be done separately for each 'study' specified in the *proj_config_template.yml*
  

## Usage notes
### Specified within nextflow_template.config:
 * Annotation files generated from VCF files will be matched by wildcard character if specified or assumes a 1-1 matching
  * e.g., if annotation file has name: **Sites_only_VCF_chr\*.vcf** and the regenie input file has name **Another_file_chr\*.pgen**, then will match based on whatever is specified by the character in the '*' position
  * e.g., if annotation file has name: **Sites_only_allchr.vcf** and the regenie input is **Another_file_allchr.bgen**, then it will generate only 1 annotation file and assumes it matches to **Another_file_allchr.bgen**
  * the wildcard character can stand-in for 1 or more alphanumeric symbols.
### For nextflow_template.config and proj_config_template.yml:
 * Any flags and values meant for Regenie (i.e., *step2_exwas_genetic* in nextflow_template.config and anything in *s1_params* and *s2_params* from the proj_config_template.yml) will be passed directly to Regenie so the flag names and the values have to be what it expects based on [Regenie documentation](https://rgcgithub.github.io/regenie/options/) (e.g., pgen/bfile have no extensions)

## Configuration files
  * exwas_pipeline.yml: conda environment file to execute the python scripts
  * nextflow_template.config: nextflow configuration
  * proj_config_template.yml: ExWAS configuration yaml files:
## program requirements (paths to be specified in proj_config_template.yml):
  * nextflow >= 23.10.0
  * python 3.10.9
  * Other required programs (plink, tabix, etc) are listed in proj_config_template.yml

## stuff to figure out...Ordered from most to least important

**Modify the ExWAS part so finds the proper annotatin ofile based on wildcard.**

**Add validation checks for configuration format with informative errors**

**simplify config file for specifying masks**

**Figure out how to put the conda environment in a docker/singularity and run with those image.**
  * how to run one container (the vep container) inside another container (the workflow container)????

**Parameters are specified in 2 files right now. Can it be merged?**
  * nextflow_template.config:
    * specifies how to run nextflow in conda environment without conda, and potentially with docker, etc
    * specifies the location of the vcf files for annotation and ExWAS.
      * specified here because using the wildcard, nextflow will pass each vcf through the pipeline in parallel (e.g., we can annotate vcf for each chr in parallel or run Regenie for each chr in parallel)
  * proj_config_template.yml:
    * specifies all other parameters to run the scripts.
      * specified here because there are certain parameters that the script expects to be like a dictionary or list, etc. Cannot be done (or not sure how to get it done) using the commandline argument.
      * Right now, it reads this yml file will store to get those data structures

**Different 'study' with different masks and annotations are specified in proj_config_template.yml so will be done sequentially.**
* e.g., a "Regeneron" and "Genomics England" ExWAS can be specified in there, but will just do one after another
  * Perhaps specify this in the nextflow_template.config so each study is done in parallel?
  * e.g., input: Regeneron chr1, ... Regeneron chr23, Genomics England chr1,...Genomics England chr23 and run these in parallel?

**Put the config files (filled in version) here?**

**Find some small public example data to test installation and stuff and write test cases.**

**Can add LDSC as a workflow**
  * Provided the user download all the files and give all the flag, should be pretty easy.

**Add functions to handle outputs from more plugins**
 * done in:     python_helpers/vep_helpers/parse_vep.py
