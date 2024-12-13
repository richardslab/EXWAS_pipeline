# number of cases far exceeds what was reported on GWAS catalog
# case/control definitions are different for cataract which showed no correlation.
# case control definition more similar with Hypothyroidisms

import os,shutil,yaml,csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

backman_cataract_gwas = pd.read_csv(
  "/home/richards/kevin.liang2/scratch/exwas_pipeline/data/backman_exwas_results/GCST90083419_buildGRCh38.tsv.gz",
  sep="\t",quoting=csv.QUOTE_NONE
)
alpha_cataract_file = pd.read_csv(
  "/scratch/richards/yiheng.chen/project14_ExWAS_AlphaMissense/data/UKB_eur_binary_phenotypes.txt",
  sep="\t",quoting=csv.QUOTE_NONE
)
alpha_cataract_file.Cataract.sum()
len(alpha_cataract_file) - alpha_cataract_file.Cataract.sum()

cataract_case = alpha_cataract_file.query("Cataract == 1")
cataract_control = alpha_cataract_file.query("Cataract == 0")

hypothyroidism_case = alpha_cataract_file.query("Hypothyroidism == 1")
hypothyroidism_control = alpha_cataract_file.query("Hypothyroidism == 0")