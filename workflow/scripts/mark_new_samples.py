#!/usr/bin/env python3

__author__ = "Boas van der Putten"
__version__ = "0.0.1"

import pandas as pd


def main(args):
  list_samples = list(args.samples)

  df = pd.read_csv(args.input, sep='\t', header=None, names=['query', 'ref', 'SNPs'])
  list_cluster = list(df['query']) + list(df['ref'])
  set_cluster = set(list_cluster)
  
  dict_marked = {}
  for strain in set_cluster:
    if strain in list_samples:
      dict_marked[strain] = 'New'
    else:
      dict_marked[strain] = 'Old'

  out_df = pd.DataFrame.from_dict(dict_marked, orient='index')
  out_df = out_df.reset_index()
  out_df.columns = ['id', 'batch']
  out_df.to_csv(args.output, index=False)


if __name__ == "__main__":
  import argparse
  
  parser = argparse.ArgumentParser(description='Convert SNP table to SNP network in dot format.')

  parser.add_argument("-i", "--input", dest="input", help="Input table (tsv format)", metavar="INPUT FILE", required=True, type=str)
  parser.add_argument("-o", "--output", dest="output", help="Output table (csv format)", metavar="OUTPUT FILE", required=True, type=str)
  parser.add_argument('samples', nargs='+', help="List of sample names of the current new batch", type=str)
  parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))

  args = parser.parse_args()
  
  main(args)
