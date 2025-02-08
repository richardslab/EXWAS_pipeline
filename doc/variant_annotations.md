## VEP Annotation notes
  * VCF for generating annotation files are specified separately from the input to run Regenie. 
      * This way, sites-only VCF can be used to generate annotation, as it is smaller and faster to run.
      * As long as the variants are the same, should not need to regenerate the annotations for each ExWAS (e.g., same set of annotations for all males and females, stratified analyses, etc) (I think...)
    * If VCF contains genetic data, a sites-only VCF will be created by the --drop-genotypes flag in bcftools. If the VCF is a sites-only VCF, then this flag will simply not have any effect (I think...)
  * Using VEP
    * The VEP image have no plugins and none of the cache files required to run any plugins. It only has vep installed.
      * Have to download everything, then specify these location in *proj_config_template.yml*
  * The plugins that are parsed right now:
     * IMPACT
     * LoFtee
     * REVEL
     * dbNSFP
       * alphamissense_pred
       * CADD_phred
       * EVE_Class25_pred
       * LRT_pred
       * MutationTaster_pred
       * Polyphen2_HDIV_pred
       * Polyphen2_HVAR_pred
       * SIFT4G_pred
       * SIFT_pred