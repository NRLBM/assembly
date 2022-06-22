#!/bin/env python3

import pandas as pd
import networkx as nx
import subprocess
import os

def clean_filenames(df, list_cols):
  for col in list_cols:
    df[col] = df[col].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
  return df

def filter_mash_dist(df, isolate, threshold):
  filename = ''.join([isolate, '.tsv'])
  filepath = '/'.join(['tmp_data', args.timestamp, 'mash_out', filename])
  tmp_df = pd.read_csv(filepath, sep='\t', names=['query', 'ref', 'distance', 'p-value', 'hashes'])
  tmp_df = clean_filenames(tmp_df, ['query', 'ref'])
  tmp_df = tmp_df[tmp_df['distance'] <= threshold]
  df = pd.concat([df, tmp_df])
  return df

def remove_graphs_without_isolates(graph, list_isolates):
  combined_graph = nx.Graph()
  for nodelist in nx.connected_components(graph):
    count_subgraph = 0
    for isolate in list_isolates:
      if isolate in nodelist:
        count_subgraph += 1
    if count_subgraph > 0:
      subgraph = graph.subgraph(nodelist)
      combined_graph = nx.compose(combined_graph, subgraph)
  return combined_graph

def filter_stricter_threshold(graph, nodelist, threshold):
  subgraph = nx.Graph(graph.subgraph(nodelist))
  subgraph.remove_edges_from([(n1, n2) for n1, n2, distance in subgraph.edges(data='distance') if distance > threshold])
  subgraph = remove_graphs_without_isolates(subgraph, nodelist)
  return subgraph

def create_subgraphs(list_isolates, threshold, max_nr_isolates):
  print(f"Setting current threshold to {threshold}")
  current_threshold = threshold
  full_df = pd.DataFrame()
  for isolate in list_isolates:
    full_df = filter_mash_dist(full_df, isolate, current_threshold)
  G = nx.from_pandas_edgelist(full_df, 'query', 'ref', 'distance')
  ## Remove graphs without isolates
  G = remove_graphs_without_isolates(G, list_isolates)
  ## Filter if subgraphs are too large
  total_graph = nx.Graph()
  for nodelist in nx.connected_components(G):
    nr_isolates = len(nodelist)
    subgraph = G.subgraph(nodelist)
    while nr_isolates > max_nr_isolates:
      current_threshold = current_threshold * 0.9
      subgraph = filter_stricter_threshold(subgraph, nodelist, current_threshold)
      subgraph = remove_graphs_without_isolates(subgraph, list_isolates)
      nr_isolates = len(subgraph.nodes())
    total_graph = nx.compose(total_graph, subgraph)
  return total_graph

def main(args):
  list_isolates = list(args.samples)

  clusters = create_subgraphs(list_isolates, args.threshold, args.max_nr_strains_cluster)
  cluster_number = 1
  if not os.path.isdir(args.out):
    os.makedirs(args.out)
  for subcluster in nx.connected_components(clusters):
    outpath = args.out + '/cluster_' + str(cluster_number) + '.txt'
    with open(outpath, 'w') as file:
      for strain in subcluster:
        assembly_name = ''.join([strain, '.fasta'])
        command = ' '.join(['find', './output/', '-type f', '-name', assembly_name, '-print'])
        subprocess_out = subprocess.run(command, capture_output=True, shell=True)
        strain_path = subprocess_out.stdout.decode("utf-8").rstrip('\n').split('\n')
        file.write('\t'.join([strain, strain_path[0], '\n']))
    cluster_number += 1

if __name__ == '__main__':
  import argparse
  
  # Parse arguments
  parser = argparse.ArgumentParser(description='Combine mash distances to get lists of strains to call SNPs for')

  parser.add_argument('samples', nargs='+', help="List of sample names to combine", type=str)
  parser.add_argument('-t', '--threshold', help="Mash threshold to filter", dest='threshold', default=0.005, type=float)
  parser.add_argument('--max-strains-cluster', help="Maximum number of strains per cluster", dest='max_nr_strains_cluster', default=100)
  parser.add_argument('-o', '--out', help="Output directory", dest='out', required=True)
  parser.add_argument('--timestamp', help="Timestamp", dest='timestamp', required=True)

  args = parser.parse_args()

  main(args)
