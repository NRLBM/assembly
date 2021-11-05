#!/usr/bin/env python3

import argparse
import json
import os
import pandas as pd

parser = argparse.ArgumentParser(description='Convert PubMLST typing (traditional MLST scheme) to table.')

parser.add_argument('--input', dest='input', required=True, type=str, help='Input JSON file with MLST allele calls from PubMLST API')
parser.add_argument('--output', dest='output', required=True, type=str, help='Output TSV file')
parser.add_argument('--loci', dest='loci', required=True, type=str, help='MLST scheme loci to use')

args = parser.parse_args()

with open(args.input) as file:
  results_dict = json.load(file)

mlst_dict = {}

isolate_name = os.path.splitext(os.path.basename(args.input))[0]
mlst_dict['Isolate'] = isolate_name

if 'fields' in results_dict:
  if 'ST' in results_dict['fields']:
    mlst_dict['ST'] = results_dict['fields']['ST']
  else:
    mlst_dict['ST'] = '-'

  if 'clonal_complex' in results_dict['fields']:
    mlst_dict['CC'] = results_dict['fields']['clonal_complex']
  else:
    mlst_dict['CC'] = '-'
else:
  mlst_dict['ST'] = '-'
  mlst_dict['CC'] = '-'

mlst_loci = args.loci.split(',')

multiple_allele_loci = []

for locus in mlst_loci:
  if locus in results_dict['exact_matches']:
    if len(results_dict['exact_matches'][locus]) > 1:
      multiple_allele_loci.append(locus)
    mlst_dict[locus] = results_dict['exact_matches'][locus][0]['allele_id']
  else:
    mlst_dict[locus] = '-'

if len(multiple_allele_loci) == 0:
  mlst_dict['multiple_alleles'] = '-'
else:
  mlst_dict['multiple_alleles'] = '|'.join(multiple_allele_loci)

df = pd.DataFrame(mlst_dict, index=[0])

df.to_csv(args.output, index=False, sep = '\t')
