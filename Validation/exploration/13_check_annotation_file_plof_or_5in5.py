"""
  Check to make sre all variants annotated as deleterious 5in5 were not IMPACT HIGH

  Based on the regeneron annotations/masks, variants are first pLoF and if not then see if they satisfy deleterious criteria

  # Results:
  For each gene, each variants, those annotated as deleterious 5 in 5 did not have a IMPACT=HIGH entry in the annotation file which makes sense as these would be considered as pLoF and pLoF takes precedence.
"""
import os,shutil,yaml,re,sqlite3,csv
import pandas as pd
import numpy as np


chr1_annotation = pd.read_csv(
  "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/plof_or_5in5/regeneron/annotations_wes_qc_chr1_sitesonly.txt",
  header=None,sep="\t"
).rename(
  {
    0:"SNP",1:"GENE",2:"annotation"
  },axis=1
)
deleterious_var = chr1_annotation.query("annotation == 'deleterious_5_of_5'")
plof_var = chr1_annotation.query("annotation == 'pLoF'")

chr1_vep_output = pd.read_csv(
  "/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/plof_or_5in5/3_annotation_results_wes_qc_chr1_sitesonly.txt.gz",
  comment="#",sep="\t",header=None,quoting=csv.QUOTE_NONE
).rename(
  {
    0:"Uploaded_variation",
    1:"Location",
    2:"Allele",
    3:"Gene",
    4:"Feature",
    5:"Feature_type",
    6:"Consequence",
    7:"cDNA_position",
    8:"CDS_position",
    9:"Protein_position",
    10:"Amino_acids",
    11:"Codons",
    12:"Existing_variation",
    13:"Extra"
  },axis=1
)


del_var_annotations = pd.merge(
  deleterious_var,
  chr1_vep_output,
  'inner',
  left_on=['SNP','GENE'],right_on=['Uploaded_variation','Gene']
)
del_var_annotations['impact'] = del_var_annotations['Extra'].apply(
  lambda val: val.split(";")[0].split("=")[1]
)
del_var_annotations.impact.unique()
del_var_annotations.query("impact == 'HIGH'")