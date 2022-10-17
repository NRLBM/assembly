#!/usr/bin/env python3

__author__ = "Boas van der Putten"
__version__ = "0.0.1"

import argparse

parser = argparse.ArgumentParser(description='Summarise assembled isolates')

parser.add_argument("--species", required=True, help="Species to select samples for", type=str)
parser.add_argument('--timestamp', required=True, help='Timestamp to parse results for', type=str)
parser.add_argument("--thresholds", dest='thresholds', help='Quality control thresholds for sample selection', type=str, default='workflow/config/qc_thresholds.json')
parser.add_argument("--exclude", dest='exclude', help='Specific samples to exclude regardless of QC', nargs='+', default = [])
parser.add_argument('--output', required=True, dest='output', help='Destination of output', type=str)
parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))

args = parser.parse_args()

import pandas as pd
import numpy as np
import json

def select_all_samples_for_species(timestamp, species):
  summary_path = '/'.join(['backup', timestamp, 'summary.csv'])
  df = pd.read_csv(summary_path)
  df_species = df.query('`Top species` == @species')
  return df_species

def read_qc_filter_thresholds(species, thresholds_file):
  with open(thresholds_file) as file:
    all_qc_thresholds = json.load(file)
  qc_thresholds_species = all_qc_thresholds[species]
  return qc_thresholds_species

def qc_filter_species(df_species, thresholds_dict):
  min_N50, min_genome_size, max_genome_size, min_sequencing_depth, max_number_contigs = thresholds_dict['min_N50'], thresholds_dict['min_genome_size'], thresholds_dict['max_genome_size'], thresholds_dict['min_sequencing_depth'], thresholds_dict['max_number_contigs']
  df_out = pd.DataFrame()
  df_out['Isolate'] = df_species['Isolate']
  df_out['N50'] = np.select([df_species['N50'] >= min_N50], ['PASS'], default='FAIL')
  df_out['Total assembly size'] = np.select([(df_species['Total assembly size'] >= min_genome_size) & (df_species['Total assembly size'] <= max_genome_size)], ['PASS'], default='FAIL')
  df_out['Sequencing depth'] = np.select([df_species['Sequencing depth'] >= min_sequencing_depth], ['PASS'], default='FAIL')
  df_out['Number of contigs'] = np.select([df_species['Number of contigs'] <= max_number_contigs], ['PASS'], default='FAIL')
  df_out['QC_decision'] = np.select([(df_out.iloc[:,1:] == 'PASS').all(axis=1)], ['PASS'], default='FAIL')
  return df_out

def remove_specific_samples(df_qc_filtered, exclusion_list):
  if len(exclusion_list) == 0:
    return df_qc_filtered
  else:
    df_qc_filtered.loc[df_qc_filtered['Sample'].isin(exclusion_list), 'QC_decision'] = 'MANUALLY_EXCLUDED'
    return df_qc_filtered

def save_output(df_out, output):
  df_out.to_csv(output, sep = '\t', index=False)

def main(args):
  df_species = select_all_samples_for_species(args.timestamp, args.species)
  thresholds_dict = read_qc_filter_thresholds(args.species, args.thresholds)
  df_qc_filtered = qc_filter_species(df_species, thresholds_dict)
  df_qc_filtered_excluded = remove_specific_samples(df_qc_filtered, args.exclude)
  save_output(df_qc_filtered_excluded, args.output)

if __name__ == '__main__':
  main(args)
