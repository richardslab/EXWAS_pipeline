## Configuration files
  * exwas_pipeline.yml: conda environment file to execute the python scripts
  * nextflow_template.config: nextflow configuration
  * proj_config_template.yml: ExWAS configuration yaml files:

###  nextflow_template.config: Specify genetic file inputs
* Variant annotation takes in VCF files and Burden testing requires PGEN, BGEN, or plink binary files (.bed/.bim/.fam)
  * Variants are expected to have the same IDs between the two sets of files
  * Matching is done by wildcard characters or expects 1-1 matching
  
  
### proj_config_template.yml: specification of study masks
  * Masks are defined in the "mask_definitions" flag in the form of a python dictionary
    ```
    "study_name":{"mask_name":["variant_annotation1","variant_annotation2",...,"variant_annotation3"]}
    ```
     * e.g.,:
        * to specify a study with 2 masks. mask1 includes all pLoF and deleterious variants and mask2 includes only pLoF:
          ```
          "studyA":{
           "mask1":["pLoF","deleterious"],
           "mask2":["pLoF"]
          }
          ```
  * How variants are annotated based on plugins are defined in the "annotation_definitions" flag in the form of a python dictionary
    ```
    "study_name":{
      "annotation1":{
        "type":{
         "plugin1":['plugin1 criteria1','plugin1 criteria2',..,'plugin X criteria X'],
         'var_consequence':['var_consequence1','var_consequence2',...,'var_consequenceX']
        }
      },
      "annotation2":{...}
    }
    ```
    * 'var_conseqence' is optional. It define what type of variants are considered. For instance, only variants tagged with "missense_variant" by VEP or "missense_variant,3UTR_region" are included. If this flag is omitted, all variants are considered.
    * 'type': can be either "all", "any", or any numeric value. This specify criteria from how many plugins are required.
    * For instance. To specify studyA where:
        * pLoF = VEP HIGH impact
        * deleterious_5in5 based on LRT, MutationTaster, Polyphen2_HDIV, Polyphen2_HDVAR, and SIFT. Only if it is annotated as deleterious by all 5 programs. Only missense variants considered. One would do
        * deleterious_1_or_more based on LRT, MutationTaster, Polyphen2_HDIV, Polyphen2_HDVAR, and SIFT. Only if it is annotated as deleterious by all 5 programs. Only missense variants considered. One would do
      ```
      "studyA":{
        "pLoF":{
          "all":{
            "IMPACT":["HIGH"]
           }
         },
        "deleterious_5in5":{
          "all":{
            "LRT_Pred":["D"],
            "MutationTaster_pred":["A","D"],
            "Polyphen2_HDIV":["D"],
            "Polyphen2_HVAR":["D"],
            "SIFT_pred":["D"]
          },
          "var_consequence":["missense_variant"]
        },
      "deleterious_1_or_more":{
          1:{
            "LRT_Pred":["D"],
            "MutationTaster_pred":["A"],
            "Polyphen2_HDIV":["D"],
            "Polyphen2_HVAR":["D"],
            "SIFT_pred":["D"]
          },
          "var_consequence":["missense_variant"]
        }
      ```
  * The order of annotation to pick for variants that satisfy multiple annotation criteria is defined in "annotation_order" flag as a form of python dictionary
    * Notes: Even if there is only 1 annotation, this is required
    * for instance, to always select "pLoF" before "deleterious_5in5" and "deleterious_5in5" before "deleterious_1_or_more"
      ```
      {"studyA":["pLoF","deleterious_5in5","deleterious_1_or_more"]}
      ```
### proj_config_template.yml: Specify Regenie parameters
 * Any flags and values meant for Regenie (i.e., anything in *s1_params* and *s2_params* from the proj_config_template.yml) will be passed directly to Regenie so the flag names and the values have to be what it expects based on [Regenie documentation](https://rgcgithub.github.io/regenie/options/)