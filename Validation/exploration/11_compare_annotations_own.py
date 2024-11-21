# the same results with the new SQL query
import os,shutil,yaml,re,glob,csv
import pandas as pd
import numpy as np


first_annotation = pd.read_csv("/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/plof_or_5in5/regeneron/old_annotation/annotations_wes_qc_chr21_sitesonly.txt",sep="\t",quoting=csv.QUOTE_NONE,header=None).rename(
  {
    0:'snp_f',1:'gene_f',2:'annotation_f'
  },axis=1
)
new_annotation = pd.read_csv("/home/richards/kevin.liang2/scratch/exwas_pipeline/results/Validation_regeneron/plof_or_5in5/regeneron/annotations_wes_qc_chr21_sitesonly.txt",sep="\t",quoting=csv.QUOTE_NONE,header=None).rename(
  {
    0:'snp_o',1:'gene_o',2:'annotation_o'
  },axis=1
)

shared = pd.merge(
  first_annotation,new_annotation,'inner',left_on=['snp_f','gene_f'],right_on=['snp_o','gene_o']
)
all(
  shared['annotation_o'] == shared['annotation_f']
)
