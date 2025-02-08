process create_mask_files {
  storeDir params.outdir
  input: 
    path annotate_vcf_logs
    tuple val(annotation_vcf), val(ofile_suffix)
    val config_file
    val wdir
  
  output:
    // wildcard for study names
    // figure out if it can specify numbers
    path "*/masks_${ofile_suffix}.txt" 
    path "4_create_masks_${ofile_suffix}.logs", emit: "log"
  
  script:
  """
  set -o pipefail
  python -u ${baseDir}/modules/Regenie_input_preparation/03_create_mask_files.py -c ${config_file} -i ${annotation_vcf} --wdir ${wdir} | tee "4_create_masks_${ofile_suffix}.logs"
  """
}