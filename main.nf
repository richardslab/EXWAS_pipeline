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
  output:
    stdout

  script:
  // baseDir references the path containing the current nextflow script
  // to obtain parent directory, use dirname of bash
  // escape characters accordingly.
    """
    #!/bin/bash
    echo \$(dirname ${baseDir})
    """
}

workflow {
  test_python | view
}