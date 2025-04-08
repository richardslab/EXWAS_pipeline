process create_mask_files {
  publishDir("${params.outdir}/REGENIE_INPUTS",mode:"link")
  errorStrategy "terminate"
  label "light","instant"

  input: 
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val anno_res
  
  output:
    // wildcard for study names
    path "**/masks_${ofile_suffix}.txt", emit: "mask_files"
    path "4_create_masks_${ofile_suffix}.logs", emit: "log"
  
  script:
  """
  set -o pipefail
  create_mask_files.py -c ${config_file} -i ${annotation_vcf} --input_dir "${params.outdir}/ANNOTATE_VARIANTS/03_variant_annotations" | tee "4_create_masks_${ofile_suffix}.logs"
  """
}