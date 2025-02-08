# EXWAS_pipeline

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## This documentations is split into 2 sections:

- For detailed usage descriptions
  - [Input/Output files](doc/input_outputs.md)
  - [Regenie Burden testing information](doc/burden_testing.md)
  - [VEP variant annotations](doc/variant_annotations.md)
  - [How to specify the coniguration](doc/configuration_file_descriptions.md)
  - [Citations](./CITATIONS.md)


## Getting started
```
nextflow run \
  ./workflows/exwas.nf \ 
  -c ./nextflow_template.config \
  -profile conda
```
> [!WARNING]
> Please make sure information in these configuration files are correct
> - [nextflow_template.config](./nextflow_template.config)
> - [proj_config_template.yml](./proj_config_template.yml)

## Pipeline descriptions:

This is a set of python/nextflow scripts that is used to run Gene-based burden testing using Regenie. This pipeline is divided into 2 parts.

**Part 1:**
This is to annotated all the variants using VEP with plugins specified by the users and generate input files to run Regenie. This includes annotation files, mask files, and the setlist files required.

The apptainer VEP image can be created using the definition files provided. Default is version 105. 

**Part 2:**
This runs Regenie step 1 and step 2 with user defined parameters. Step 1 expects 1 set of plink files of genotyped variants. This step is performed once. Step 2 is per analysis (as specified by the user in the configuration files)

## program requirements (paths to be specified in proj_config_template.yml):
  * nextflow >= 23.10.0
  * python 3.10.9
  * SQLITE3 3.40.1
  * Other required programs (plink, tabix, etc) are listed in proj_config_template.yml
  * Python packages required will be specified via the CONDA environment

## General Notes
  * Making conda environment on first run will take time. As long as the conda cache is not deleted, the environment will not be made again.
  * For all steps in the pipeline, as long as the log files are not deleted, the steps will be skipped.
    * **Notes** this also means that if the log files are moved from default location or renamed, the steps will be re-ran and files will be overwritten.
  * Remember to set this flag to *bcftools_param_set_id* 0 to use the original IDs in the VCF file.
    * if set to 1, then the annotation files will be using the modified SNP id that is chr:pos:ref:alt **after left alignment**. Needs to take extra step to ensure this modified ID matches your ExWAS input data if using this flag.


## Citations
If you use Exwas_pipeline for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

