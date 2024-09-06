"""Parse the consequences of the vep file for the plugins that we used
"""

# https://usf.app.box.com/s/p2etici4mp5noboju6g0mnqfjx8ymwzd

# alphamissense order
## P: Pathogenic
## B: Benign
## A: Ambiguous
ALPHAMISSENSE_ORDER = ('P',"B","A")

# EVE class 25 order
## P: pathogenic
## B: benign
## U: uncertain
EVE_CLASS25_ORDER = ("P","B","U")

# LRT_pred order
## D: deleterious
## N: neutral
## U: unknown
LRT_ORDER = ("D","N","U")

# mutation taster order
#: automatic before non-automatic https://www.mutationtaster.org/info/documentation.html
## A: disease_causing_automatic
## D: Disease causing
## N: Polymorphism
## P: Polymorishm automatic
MUTATIONTASTER_ORDER = ("A","D","P","N")


def __parse_IMPACT(consequence,consequence_elem):
   val = consequence_elem[1]
   assert(len(consequence) == 1),f"invalid IMPACT results {consequence}"
   return val

def __parse_alphamissense(consequence,consequence_elem):
  alpha_preds = [x.upper().strip() for x in consequence_elem[1].split(",")]
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

def __parse_eve_class25(consequence,consequence_elem):
  eve_class25_preds = [x.upper().strip() for x in consequence_elem[1]]
  assert(
    all(
      [x in EVE_CLASS25_ORDER for x in eve_class25_preds]
    )
  ),f"invalid EVE class25 pred {consequence}"
  for i in EVE_CLASS25_ORDER:
    if i in eve_class25_preds:
      return i
  assert(False)

def __parse_mutationtaster(consequence,consequence_elem):
  mutationtaster_preds = [x.upper().strip() for x in consequence_elem[1]]
  assert(
    all(
      [x in MUTATIONTASTER_ORDER for x in mutationtaster_preds]
    )
  ),f"invalid EVE class25 pred {consequence}"
  for i in MUTATIONTASTER_ORDER:
    if i in mutationtaster_preds:
      return i
  assert(False)

def parse_var_consequence(annotation):
  """summarize the consequences for each variant for the plugins that we have

  Args:
      annotation (str): 1 entry of the 'extra' columns
  """

  var_consequences = dict()
  all_consequences = annotation.split(";")
  for consequence in all_consequences:
    consequence_elem = [x.strip() for x in consequence.split("=")]
    if consequence_elem[0] == "IMPACT":
      var_consequences['IMPACT']=__parse_IMPACT(consequence,consequence_elem)

    elif consequence_elem[0] == 'AlphaMissense_pred':
      var_consequences['AlphaMissense_pred']=__parse_alphamissense(consequence,consequence_elem)

    elif consequence_elem[0] == 'CADD_phred':
      var_consequences['CADD_phred'] = __parse_CADD(consequence_elem)

    elif consequence_elem[0] == 'EVE_Class25_pred':
      var_consequences['EVE_Class25_pred'] = __parse_eve_class25(consequence,consequence_elem)
    
    elif consequence_elem[0] == 'MutationTaster_pred':
      var_consequences['MutationTaster_pred'] = __parse_mutationtaster(consequence,consequence_elem)
    
  return