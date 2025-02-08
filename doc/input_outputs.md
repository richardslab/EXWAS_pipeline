## Usage notes
### Input
**Common inputs**
 * project configuration yaml file
 * nextflow configuration file
   
**Part 1**
 * 1 or more VCF files for annotation purposes

**Part 2**
 * Part 1 outputs
 * 1 or more genetic dataset for Regenie. Format can be pgen, bgen, or bed (same as what Regenie takes)

### Output
**Part 1**
 * sites only VCF for each VCF files provided
 * VEP annotation output and summary information
 * Summarized output results in an indexed SQL database
 * Mask files, annotation files, and setlist files for each analysis specified by the user
   
**Part 2**
 * Regenie S1 output
 * Regenie S2 output