#!/usr/bin/env python3

import json
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Predict vaccine coverage based on gMATS method')

parser.add_argument('--input', dest='input', required=True, type=str, help='Input JSON file from PubMLST typing')
parser.add_argument('--output', dest='output', required=True, type=str, help='TSV output file')

args = parser.parse_args()

isolate_name = os.path.splitext(os.path.basename(args.input))[0]

with open(args.input) as file:
  pubmlst_typing = json.load(file)['exact_matches']

tuple_fHbp_peptide_covered = ('1', '2', '4', '14', '15', '37', '89', '90', '110', '144', '224', '232', '245', '249', '252', '510')
tuple_fHbp_peptide_uncovered = ('16', '19', '21', '22', '23', '24', '25', '29', '30', '31', '45', '47', '106', '109', '119', '160', '174', '213')

tuple_NHBA_peptide_covered = ('1', '2', '3', '5', '10', '20', '21', '113', '243')
tuple_NHBA_peptide_uncovered = ('6', '13', '17', '18', '19', '24', '25', '30', '31', '43', '47', '58', '112', '114', '120', '122', '160', '187', '253')

def check_result(typing_result, peptide, tmp_tuple):
  result = False
  typing_list = typing_result[peptide]
  result = sum([x['allele_id'] in tmp_tuple for x in typing_list])
  alleles = '|'.join([x['allele_id'] for x in typing_list])
  return result, alleles

def make_decision(covered, uncovered):
  if covered:
    return 'covered'
  elif uncovered:
    return 'uncovered'
  else:
    return 'unpredictable'

def make_decision_porA(typing_result):
  decision = 'uncovered'
  typing_list = typing_result['PorA_VR2']
  if sum([x['allele_id'] == '4' for x in typing_list]) > 0:
    decision = 'covered'
  alleles = '|'.join([x['allele_id'] for x in typing_list])
  return decision, alleles

if 'fHbp_peptide' in pubmlst_typing:
  fHbp_covered, fHbp_alleles = check_result(pubmlst_typing, 'fHbp_peptide', tuple_fHbp_peptide_covered)
  fHbp_uncovered, fHbp_alleles = check_result(pubmlst_typing, 'fHbp_peptide', tuple_fHbp_peptide_uncovered)
  fHbp_decision = make_decision(fHbp_covered, fHbp_uncovered)
else:
  fHbp_decision = 'missing'
  fHbp_alleles = '-'

if 'NHBA_peptide' in pubmlst_typing:
  NHBA_covered, NHBA_alleles = check_result(pubmlst_typing, 'NHBA_peptide', tuple_NHBA_peptide_covered)
  NHBA_uncovered, NHBA_alleles = check_result(pubmlst_typing, 'NHBA_peptide', tuple_NHBA_peptide_uncovered)
  NHBA_decision = make_decision(NHBA_covered, NHBA_uncovered)
else:
  NHBA_decision = 'missing'
  NHBA_alleles = '-'

if 'PorA_VR2' in pubmlst_typing:
  PorA_VR2_decision, PorA_VR2_alleles = make_decision_porA(pubmlst_typing)
else:
  PorA_VR2_decision = 'missing'
  PorA_VR2_alleles = '-'

list_decisions = [fHbp_decision, NHBA_decision, PorA_VR2_decision]

if 'covered' in [fHbp_decision, NHBA_decision, PorA_VR2_decision]:
  total_coverage = 'covered'
elif ('missing' in list_decisions) or ('unpredictable' in list_decisions):
  total_coverage = 'unpredictable'
else:
  total_coverage = 'uncovered'

df = pd.DataFrame({'Isolate': isolate_name, 'Bexsero_coverage': total_coverage, 'fHbp_coverage': fHbp_decision, 'NHBA_coverage': NHBA_decision, 'PorA_VR2_coverage': PorA_VR2_decision,
                   'fHbp_allele': fHbp_alleles, 'NHBA_allele': NHBA_alleles, 'PorA_VR2_allele': PorA_VR2_alleles}, index=[0])
df.to_csv(args.output, sep = '\t', index=False)
