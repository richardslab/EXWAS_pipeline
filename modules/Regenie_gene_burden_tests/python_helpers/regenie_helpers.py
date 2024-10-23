import re,glob,os
from pathlib import Path
def __find_regenie_supplementary_files(input_vcf,nxtflow_g,nxtflow_annotation,study_dir):

  # if wildcard present extract wildcard character
  if "*" in nxtflow_g:
    wildcard_characters = re.sub(re.sub("\*","",nxtflow_g),"",input_vcf)
    assert(
      re.sub(wildcard_characters,"",input_vcf) == re.sub("\*","",nxtflow_g)
    ),f"issue extracting wildcard characters"
  else:
    wildcard_characters = None
  
  if wildcard_characters is not None:
    expected_vcf_name = re.sub("\*",wildcard_characters,Path(nxtflow_annotation).stem)
  else:
    expected_vcf_name = Path(nxtflow_annotation).stem
  
  expected_annotation_file = os.path.join(study_dir,f"annotations_{expected_vcf_name}.txt")
  expected_mask_file = os.path.join(study_dir,f"masks_{expected_vcf_name}.txt")
  expected_setlist_file = os.path.join(os.path.dirname(study_dir),f"6_{expected_vcf_name}.setlist")

  assert(
    os.path.isfile(expected_annotation_file) and
    os.path.isfile(expected_mask_file) and 
    os.path.isfile(expected_setlist_file)
  ),f"Not all expected files present\nExpected {expected_annotation_file}\n{expected_mask_file}\n{expected_setlist_file}"

  return expected_annotation_file,expected_mask_file,expected_setlist_file,wildcard_characters