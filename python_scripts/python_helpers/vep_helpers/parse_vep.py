"""Parse the consequences of the vep file for the plugins that we used

Will aggregate all annotations across transcripts picking the most deleterious consequence for each plugin.

e.g., alphamissesnse may have 6 consequence for 6 transcripts for each gene (B,B,A,P,B,B), the result will be P for this variant based on it's predicted effect on the 4th transcript.

Same is done for all other plugins.

Will also assert that the annotations are the same across all plugins
"""


def __parse_sift4g_pred(consequence,consequence_elem,SIFT_ORDER):
  sift_pred = [x.upper().strip() for x in consequence_elem[1].split(",")]
  sift_pred = [x for x in sift_pred if x!='.']
  if len(sift_pred) == 0:
    return None
  assert(
    all(
      [x in SIFT_ORDER for x in sift_pred]
    )
  ),f"invalid SIFT results {consequence}"

  # the first occurence is the most severe
  for i in SIFT_ORDER:
    if i in sift_pred:
      return i
  assert False

def __parse_sift4g_pred(consequence,consequence_elem,SIFT4G_ORDER):
  sift4g_pred = [x.upper().strip() for x in consequence_elem[1].split(",")]
  sift4g_pred = [x for x in sift4g_pred if x!='.']
  if len(sift4g_pred) == 0:
    return None
  assert(
    all(
      [x in SIFT4G_ORDER for x in sift4g_pred]
    )
  ),f"invalid SIFT4G results {consequence}"

  # the first occurence is the most severe
  for i in SIFT4G_ORDER:
    if i in sift4g_pred:
      return i
  assert False


def __parse_polyphen2_hvar(consequence,consequence_elem,POLYPHEN2_HVAR_ORDER):
  p2hvar_order = [x.upper().strip() for x in consequence_elem[1].split(",")]
  p2hvar_order = [x for x in p2hvar_order if x!='.']
  if len(p2hvar_order) == 0:
    return None
  assert(
    all(
      [x in POLYPHEN2_HVAR_ORDER for x in p2hvar_order]
    )
  ),f"invalid Polyphen2 HVAR results {consequence}"

  # the first occurence is the most severe
  for i in POLYPHEN2_HVAR_ORDER:
    if i in p2hvar_order:
      return i
  assert False


def __parse_polyphen2_hdiv(consequence,consequence_elem,POLYPHEN2_HDIV_ORDER):
  p2hdiv_order = [x.upper().strip() for x in consequence_elem[1].split(",")]
  p2hdiv_order = [x for x in p2hdiv_order if x!='.']
  if len(p2hdiv_order) == 0:
    return None
  assert(
    all(
      [x in POLYPHEN2_HDIV_ORDER for x in p2hdiv_order]
    )
  ),f"invalid Polyphen2 HDIV results {consequence}"

  # the first occurence is the most severe
  for i in POLYPHEN2_HDIV_ORDER:
    if i in p2hdiv_order:
      return i
  assert False

def __parse_lrt_order(consequence,consequence_elem,LRT_ORDER):
  lrt_preds = consequence_elem[1].upper()
  assert(
    len(consequence_elem) == 2 and
    lrt_preds in LRT_ORDER
  ),f"invalid LRT_pred results {consequence}"
  
  return lrt_preds


def __parse_lof_order(consequence,consequence_elem,LOF_ORDER):
  lof_pred = consequence_elem[1].upper()
  assert(
    len(consequence_elem)==2 and
    lof_pred in LOF_ORDER
  ),f"invalid LoF result {consequence}"
  return lof_pred

def __parse_IMPACT(consequence,consequence_elem,IMPACT_ORDER):
   val = consequence_elem[1]
   assert(
     len(consequence_elem) == 2 and 
     val in IMPACT_ORDER
    ),f"invalid IMPACT results {consequence}"
   return val

def __parse_alphamissense(consequence,consequence_elem,ALPHAMISSENSE_ORDER):
  alpha_preds = [x.upper().strip() for x in consequence_elem[1].split(",")]
  alpha_preds = [x for x in alpha_preds if x!='.']
  if len(alpha_preds) == 0:
    return None
  assert(
    all(
      [x in ALPHAMISSENSE_ORDER for x in alpha_preds]
    )
  ),f"invalid Alphamissense prediction {consequence}"
  for i in ALPHAMISSENSE_ORDER:
    if i in alpha_preds:
      return i
  assert(False)

def __parse_CADD(consequence_elem):
  return(float(consequence_elem[1]))

def __parse_eve_class25(consequence,consequence_elem,EVE_CLASS25_ORDER):
  eve_class25_preds = [x.upper().strip() for x in consequence_elem[1].split(',')]
  eve_class25_preds = [x for x in eve_class25_preds if x!='.']
  if len(eve_class25_preds) == 0:
    return None
  
  assert(
    all(
      [x in EVE_CLASS25_ORDER for x in eve_class25_preds]
    )
  ),f"invalid EVE class25 pred {consequence}"
  for i in EVE_CLASS25_ORDER:
    if i in eve_class25_preds:
      return i
  assert(False)

def __parse_mutationtaster(consequence,consequence_elem,MUTATIONTASTER_ORDER):
  mutationtaster_preds = [x.upper().strip() for x in consequence_elem[1].split(",")]
  mutationtaster_preds = [x for x in mutationtaster_preds if x!='.']
  if len(mutationtaster_preds) == 0:
    return None
  assert(
    all(
      [x in MUTATIONTASTER_ORDER for x in mutationtaster_preds]
    )
  ),f"invalid EVE class25 pred {consequence}"
  for i in MUTATIONTASTER_ORDER:
    if i in mutationtaster_preds:
      return i
  assert(False)

def parse_var_consequence(annotation,CONST):
  """summarize the consequences for each variant for the plugins that we have

  Args:
      annotation (str): 1 entry of the 'extra' columns
  """

  var_consequences = dict()
  all_consequences = annotation.split(";")
  for consequence in all_consequences:
    consequence_elem = [x.strip() for x in consequence.split("=")]
    if consequence_elem[0] == "IMPACT":
      var_consequences['IMPACT']=__parse_IMPACT(consequence,consequence_elem,CONST['IMPACT'])

    elif consequence_elem[0] == "LoF":
      var_consequences["LoF"] = __parse_lof_order(consequence,consequence_elem,CONST['LoF'])

    elif consequence_elem[0] == 'AlphaMissense_pred':
      var_consequences['AlphaMissense_pred']=__parse_alphamissense(consequence,consequence_elem,CONST['AlphaMissense_pred'])

    elif consequence_elem[0] == 'CADD_phred':
      var_consequences['CADD_phred'] = __parse_CADD(consequence_elem)

    elif consequence_elem[0] == 'EVE_Class25_pred':
      var_consequences['EVE_Class25_pred'] = __parse_eve_class25(consequence,consequence_elem,CONST['EVE_Class25_pred'])
    
    elif consequence_elem[0] == 'LRT_pred':
      var_consequences['LRT_pred'] = __parse_lrt_order(consequence,consequence_elem,CONST['LRT_pred'])
    
    elif consequence_elem[0] == 'MutationTaster_pred':
      var_consequences['MutationTaster_pred'] = __parse_mutationtaster(consequence,consequence_elem,CONST['MutationTaster_pred'])

    elif consequence_elem[0] == 'Polyphen2_HDIV_pred':
      var_consequences['Polyphen2_HDIV_pred'] = __parse_polyphen2_hdiv(consequence,consequence_elem,CONST['Polyphen2_HDIV_pred'])
    
    elif consequence_elem[0] == 'Polyphen2_HVAR_pred':
      var_consequences['Polyphen2_HVAR_pred'] = __parse_polyphen2_hvar(consequence,consequence_elem,CONST['Polyphen2_HVAR_pred'])

    elif consequence_elem[0] == 'SIFT4G_pred':
      var_consequences['SIFT4G_pred'] = __parse_sift4g_pred(consequence,consequence_elem,CONST['SIFT4G_pred'])
    
    elif consequence_elem[0] == 'SIFT_pred':
      var_consequences['SIFT_pred'] = __parse_sift4g_pred(consequence,consequence_elem,CONST['SIFT_pred'])
    
  return var_consequences