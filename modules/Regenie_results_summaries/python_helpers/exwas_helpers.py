import numpy as np
from scipy import stats
import gzip,os,yaml

def __get_study_exwas_paths(wdir,each_study):
  """helper function to get the dictionary of all the exwas results for each study.
  
  Since this is done repeatedly, best to keep it separate.

  Args:
      wdir (str): work directory (base for where eveything is)
      each_study (str): the study name as defined on proj_config

  Returns:
      dict,str: the dictionary and path to output
  """
  study_dir = os.path.join(wdir,each_study)
  regenie_s2_dir = os.path.join(study_dir,'Regenie_S2')
  summary_out = os.path.join(study_dir,'Regenie_Summaries')
  summary_file = os.path.join(summary_out,'Phenotypes_results_paths.yaml.gz')
  assert(os.path.isfile(summary_file)),f"Missing phenotype paths. Run 00_find_data.py or error in earlier steps"
  with gzip.open(summary_file,'rt') as ptr:
    study_regenie_result_paths = yaml.safe_load(ptr)
  return study_regenie_result_paths,summary_out

def compute_lambda(p_values):
  """Computes the genomic inflation factor (also known as genomic controls, or lambda) from GWAS or ExWAS p-values

  as per 
    * https://en.wikipedia.org/wiki/Genomic_control
    * http://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html

  Formula:
    lambda = median(observed chisq statistics)/(median of chi-square distribution with 1 degrees of freedom)
    

  Args:
      p_values (numpy array): an array of p_values
  """

  p_values = np.asarray(p_values)
  p_values = p_values[~np.isnan(p_values)]
  assert(len(p_values) > 0),f"No valid p-values found"

  # obtain the chisq statistics with a 50% probability under the chi-square distribution
  median_expected_chisq = stats.chi2.ppf(0.5,df=1)
  
  # for each p-value, obtain the chi-square statistics
  observed_chisq = [stats.chi2.cdf(x,df=1) for x in p_values]

  median_observed_chisq = np.median(observed_chisq)
  lambdaa = median_observed_chisq/median_expected_chisq

  return lambdaa