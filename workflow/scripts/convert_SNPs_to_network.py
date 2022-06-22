#!/usr/bin/env python3

import pandas as pd
import networkx as nx
import pydot

def main(args):
  df = pd.read_csv(args.input, sep="\t", header=None, names=['query', 'ref', 'SNPs'])
  df = df[df['query'] != df['ref']]
  G = nx.from_pandas_edgelist(df, 'query', 'ref', 'SNPs')
  G.remove_edges_from([(n1, n2) for n1, n2, distance in G.edges(data='SNPs') if distance > args.threshold])
  nx.drawing.nx_pydot.write_dot(G, args.output)

if __name__ == '__main__':
  import argparse
  
  parser = argparse.ArgumentParser(description='Convert SNP table to SNP network in dot format.')
  
  parser.add_argument("-i", "--input", dest="input", help="Input file (tsv format)", metavar="INPUT FILE", required=True, type=str)
  parser.add_argument("-o", "--output", dest="output", help="Output file (dot format)", metavar="OUTPUT FILE", required=True, type=str)
  parser.add_argument("-t", "--threshold", dest="threshold", default=25, metavar="THRESHOLD", help="SNP threshold (default: 25)", type=int)
  
  args = parser.parse_args()
  
  main(args)
