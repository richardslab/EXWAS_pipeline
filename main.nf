/*
  This is the main ExWAS workflow using Regenie.

  Consists of 2 main workflow:
    1. preparing input files for Regenie
    2. Running regenie
  
  Main data analysis are done in python.
  Data inputs are specified in a YAML parameter_file.

  workflow 1: preparing regenie input files
    #' VCF processing
    1. Align input VCF file
    2. Annotate input VCF file with VEP

    #' Regenie inputs
    3. Create mask file for the VCF
    4. Create annotation file for VCF
    5. Create set list file for VCF
*/


process test_python {
  // This specifies where the outputs are written to by the process
  storeDir "${params.outdir}/test_output"

  output:
    path "hi2.txt"

  script:
  // baseDir references the path containing the current nextflow script
    """
    python ${baseDir}/test_py/x1_test.py -c ${params.config_file}
    """
}

workflow {
  test_python | view
}