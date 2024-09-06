"""
Check the variables specified in the config file to make sure it has what it needs to run the scripts
"""

#%% General
# check all the programs exists

# check the helper files (i.e., the scripts in own package exists)

# check genome build is GRCh38 in the specification (everything uses this)

#%% annotations
# checks the that all the categories for mask_names are defined in mask_definitions
## mask_names gives tags for how the variants are annotated
## these are based on plugins that are defined in mask_definitions
## so they should all be defined
### e.g., CADD>20 means check CADD_phred to have at least 20


# check all the plugin used to define annotations are in the vep annotations


# check the vep plugin orders are in the vep annotation files

# check the numeric plugin specifications are either higher or lower
