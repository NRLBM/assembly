#!/usr/bin/env python3

import json
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Predict vaccine coverage based on gMATS method')

parser.add_argument('--input', dest='input', required=True, type=str, help='Input JSON file from PubMLST typing')
parser.add_argument('--output', dest='output', required=True, type=str, help='TSV output file')

args = parser.parse_args()


def get_scheme():
  scheme_dict = {}
  scheme_dict['gMATS_fHbp_peptide_covered'] = ('1', '2', '4', '14', '15', '37', '89', '90', '110', '144',
                                               '224', '232', '245', '249', '252', '510')
  scheme_dict['gMATS_fHbp_peptide_uncovered'] = ('16', '19', '21', '22', '23', '24', '25', '29', '30', '31',
                                                 '45', '47', '106', '109', '119', '160', '174', '213')
  scheme_dict['gMATS_NHBA_peptide_covered'] = ('1', '2', '3', '5', '10', '20', '21', '113', '243')
  scheme_dict['gMATS_NHBA_peptide_uncovered'] = ('6', '13', '17', '18', '19', '24', '25', '30', '31', '43',
                                                 '47', '58', '112', '114', '120', '122', '160', '187', '253')
  scheme_dict['BAST_fHbp_peptide_exact_match'] = ('1')
  scheme_dict['BAST_fHbp_peptide_cross_reactive'] = ('4', '10', '12', '14', '15', '37', '110', '144', '215', '232')
  scheme_dict['BAST_fHbp_peptide_none'] = ('16', '19', '21', '22', '24', '25', '29', '30', '31', '45', '47',
                                           '59', '76', '109', '119')
  scheme_dict['BAST_NadA_peptide_exact_match'] = ('8')
  scheme_dict['BAST_NadA_peptide_cross_reactive'] = ('3', '6')
  scheme_dict['BAST_NadA_peptide_none'] = ('1', '21', '100')
  scheme_dict['BAST_NHBA_peptide_exact_match'] = ('2')
  scheme_dict['BAST_NHBA_peptide_cross_reactive'] = ('1', '5', '10', '113', '243', '607')
  scheme_dict['BAST_NHBA_peptide_none'] = ('6', '9', '17', '18', '25', '30', '31', '43', '47', '63',
                                           '112', '120', '160', '187', '197')
  return scheme_dict


def process_fHbp(typing_result, results_dict, scheme_dict):
  if 'fHbp_peptide' not in typing_result:
    results_dict['fHbp_peptide', 'BAST'] = 'missing/novel allele'
    results_dict['fHbp_peptide', 'gMATS'] = 'missing/novel allele'
    results_dict['fHbp_peptide', 'alleles'] = '-'
  else:
    typing_list = typing_result['fHbp_peptide']
    covered_gMATS = sum([x['allele_id'] in scheme_dict['gMATS_fHbp_peptide_covered'] for x in typing_list]) > 0
    uncovered_gMATS = sum([x['allele_id'] in scheme_dict['gMATS_fHbp_peptide_uncovered'] for x in typing_list]) > 0
    if covered_gMATS:
      results_dict['fHbp_peptide', 'gMATS'] = 'covered'
    elif uncovered_gMATS:
      results_dict['fHbp_peptide', 'gMATS'] = 'uncovered'
    else:
      results_dict['fHbp_peptide', 'gMATS'] = 'unpredictable'
    exact_match_BAST = sum([x['allele_id'] in scheme_dict['BAST_fHbp_peptide_exact_match'] for x in typing_list]) > 0
    cross_reactive_BAST = sum([x['allele_id'] in scheme_dict['BAST_fHbp_peptide_cross_reactive'] for x in typing_list]) > 0
    none_result_BAST = sum([x['allele_id'] in scheme_dict['BAST_fHbp_peptide_none'] for x in typing_list]) > 0
    if exact_match_BAST:
      results_dict['fHbp_peptide', 'BAST'] = 'exact_match'
    elif cross_reactive_BAST:
      results_dict['fHbp_peptide', 'BAST'] = 'cross-reactive'
    elif none_result_BAST:
      results_dict['fHbp_peptide', 'BAST'] = 'none'
    else:
      results_dict['fHbp_peptide', 'BAST'] = 'insufficient data'
    results_dict['fHbp_peptide', 'alleles'] = '|'.join([x['allele_id'] for x in typing_list])
  return results_dict


def process_NHBA(typing_result, results_dict, scheme_dict):
  if 'NHBA_peptide' not in typing_result:
    results_dict['NHBA_peptide', 'BAST'] = 'missing/novel allele'
    results_dict['NHBA_peptide', 'gMATS'] = 'missing/novel allele'
    results_dict['NHBA_peptide', 'alleles'] = '-'
  else:
    typing_list = typing_result['NHBA_peptide']
    covered_gMATS = sum([x['allele_id'] in scheme_dict['gMATS_NHBA_peptide_covered'] for x in typing_list]) > 0
    uncovered_gMATS = sum([x['allele_id'] in scheme_dict['gMATS_NHBA_peptide_uncovered'] for x in typing_list]) > 0
    if covered_gMATS:
      results_dict['NHBA_peptide', 'gMATS'] = 'covered'
    elif uncovered_gMATS:
      results_dict['NHBA_peptide', 'gMATS'] = 'uncovered'
    else:
      results_dict['NHBA_peptide', 'gMATS'] = 'unpredictable'
    exact_match_BAST = sum([x['allele_id'] in scheme_dict['BAST_NHBA_peptide_exact_match'] for x in typing_list]) > 0
    cross_reactive_BAST = sum([x['allele_id'] in scheme_dict['BAST_NHBA_peptide_cross_reactive'] for x in typing_list]) > 0
    none_result_BAST = sum([x['allele_id'] in scheme_dict['BAST_NHBA_peptide_none'] for x in typing_list]) > 0
    if exact_match_BAST:
      results_dict['NHBA_peptide', 'BAST'] = 'exact_match'
    elif cross_reactive_BAST:
      results_dict['NHBA_peptide', 'BAST'] = 'cross-reactive'
    elif none_result_BAST:
      results_dict['NHBA_peptide', 'BAST'] = 'none'
    else:
      results_dict['NHBA_peptide', 'BAST'] = 'insufficient data'
    results_dict['NHBA_peptide', 'alleles'] = '|'.join([x['allele_id'] for x in typing_list])
  return results_dict


def process_NadA(typing_result, results_dict, scheme_dict):
  if 'NadA_peptide' not in typing_result:
    results_dict['NadA_peptide', 'BAST'] = 'missing/novel allele'
    results_dict['NadA_peptide', 'gMATS'] = 'uncovered'
    results_dict['NadA_peptide', 'alleles'] = '-'
  else:
    typing_list = typing_result['NadA_peptide']
    results_dict['NadA_peptide', 'gMATS'] = 'uncovered'
    exact_match_BAST = sum([x['allele_id'] in scheme_dict['BAST_NadA_peptide_exact_match'] for x in typing_list]) > 0
    cross_reactive_BAST = sum([x['allele_id'] in scheme_dict['BAST_NadA_peptide_cross_reactive'] for x in typing_list]) > 0
    none_result_BAST = sum([x['allele_id'] in scheme_dict['BAST_NadA_peptide_none'] for x in typing_list]) > 0
    if exact_match_BAST:
      results_dict['NadA_peptide', 'BAST'] = 'exact_match'
    elif cross_reactive_BAST:
      results_dict['NadA_peptide', 'BAST'] = 'cross-reactive'
    elif none_result_BAST:
      results_dict['NadA_peptide', 'BAST'] = 'none'
    else:
      results_dict['NadA_peptide', 'BAST'] = 'insufficient data'
    results_dict['NadA_peptide', 'alleles'] = '|'.join([x['allele_id'] for x in typing_list])
  return results_dict


def process_PorA_VR2(typing_result, results_dict, scheme_dict):
  if 'PorA_VR2' not in typing_result:
    results_dict['PorA_VR2', 'BAST'] = 'missing/novel allele'
    results_dict['PorA_VR2', 'gMATS'] = 'missing/novel allele'
    results_dict['PorA_VR2', 'alleles'] = '-'
  else:
    typing_list = typing_result['PorA_VR2']
    exact_match_both = sum([x['allele_id'] == '4' for x in typing_list]) > 0
    if exact_match_both:
      results_dict['PorA_VR2', 'gMATS'] = 'covered'
      results_dict['PorA_VR2', 'BAST'] = 'exact_match'
    else:
      results_dict['PorA_VR2', 'gMATS'] = 'uncovered'
      results_dict['PorA_VR2', 'BAST'] = 'none'
    results_dict['PorA_VR2', 'alleles'] = '|'.join([x['allele_id'] for x in typing_list])
  return results_dict


isolate_name = os.path.splitext(os.path.basename(args.input))[0]

with open(args.input) as file:
  pubmlst_typing = json.load(file)['exact_matches']

scheme_dict = get_scheme()

results_dict = {}
results_dict = process_fHbp(pubmlst_typing, results_dict, scheme_dict)
results_dict = process_NHBA(pubmlst_typing, results_dict, scheme_dict)
results_dict = process_NadA(pubmlst_typing, results_dict, scheme_dict)
results_dict = process_PorA_VR2(pubmlst_typing, results_dict, scheme_dict)

list_decisions_gMATS = list(filter(None, [value if key[1] == 'gMATS' else None for key, value in results_dict.items()]))
list_decisions_BAST = list(filter(None, [value if key[1] == 'BAST' else None for key, value in results_dict.items()]))

if 'covered' in list_decisions_gMATS:
  decision_gMATS = 'covered'
elif all([antigen_decision == 'uncovered' for antigen_decision in list_decisions_gMATS]):
  decision_gMATS = 'uncovered'
else:
  decision_gMATS = 'unpredictable'

if 'exact_match' in list_decisions_BAST:
  decision_BAST = 'exact_match'
elif 'cross-reactive' in list_decisions_BAST:
  decision_BAST = 'cross-reactive'
elif all([antigen_decision == 'none' for antigen_decision in list_decisions_BAST]):
  decision_BAST = 'none'
else:
  decision_BAST = 'insufficient data'

df = pd.DataFrame({'Isolate': isolate_name,
                   'Bexsero_coverage_gMATS': decision_gMATS,
                   'Bexsero_coverage_BAST': decision_BAST,
                   'fHbp_coverage_gMATS': results_dict['fHbp_peptide', 'gMATS'],
                   'NHBA_coverage_gMATS': results_dict['NHBA_peptide', 'gMATS'],
                   'NadA_coverage_gMATS': results_dict['NadA_peptide', 'gMATS'],
                   'PorA_VR2_coverage_gMATS': results_dict['PorA_VR2', 'gMATS'],
                   'fHbp_coverage_BAST': results_dict['fHbp_peptide', 'BAST'],
                   'NHBA_coverage_BAST': results_dict['NHBA_peptide', 'BAST'],
                   'NadA_coverage_BAST': results_dict['NadA_peptide', 'BAST'],
                   'PorA_VR2_coverage_BAST': results_dict['PorA_VR2', 'BAST'],
                   'fHbp_allele': results_dict['fHbp_peptide', 'alleles'],
                   'NHBA_allele': results_dict['NHBA_peptide', 'alleles'],
                   'NadA_allele': results_dict['NadA_peptide', 'alleles'],
                   'PorA_VR2_allele': results_dict['PorA_VR2', 'alleles']},
                  index=[0])


# Order dataframe columns as this is not guaranteed for Python<3.7.4
df = df[["Isolate", "Bexsero_coverage_gMATS", "Bexsero_coverage_BAST", "fHbp_coverage_gMATS",
         "NHBA_coverage_gMATS", "NadA_coverage_gMATS", "PorA_VR2_coverage_gMATS", "fHbp_coverage_BAST",
         "NHBA_coverage_BAST", "NadA_coverage_BAST", "PorA_VR2_coverage_BAST", "fHbp_allele",
         "NHBA_allele", "NadA_allele", "PorA_VR2_allele"]]

df.to_csv(args.output, sep='\t', index=False)
