#!/usr/bin/env python3

import argparse
import sys

from .fastqc import fastqc
from .species import species
from .plot import plot
from .combinetyping import combinetyping

def main():
  args = parse_args(sys.argv[1:])

  if args.subparser_name == 'species':
    species(args)

  elif args.subparser_name == 'fastqc':
    fastqc(args)

  elif args.subparser_name == 'plot':
    plot(args)

  elif args.subparser_name == 'combinetyping':
    combinetyping(args)

def parse_args(args):
  parser = argparse.ArgumentParser(description='Summarise output from WGS runs for the NRLBM pipeline')

  subparsers = parser.add_subparsers(title='Commands', dest='subparser_name')
  # def subparse for every command below this def
  species_subparser(subparsers)
  fastqc_subparser(subparsers)
  plot_subparser(subparsers)
  combinetyping_subparser(subparsers)

  if len(args) == 0:
    parser.print_help(file=sys.stderr)
    sys.exit(1)

  return parser.parse_args(args)

def species_subparser(subparsers):
  group = subparsers.add_parser('species', description='Summarise species', help='Summarise species')

  input_args = group.add_argument_group('Input arguments')
  input_args.add_argument("-i", "--input-dir", dest="inputdir", help="Location of assembly working directory and glob all summaries", metavar="HEAD DIR", default='~/assembly', type=str)
  input_args.add_argument("--selected-summaries", dest="selected_summaries", help="List selected summaries to analyse", metavar="SUMMARIES", type=str, nargs='+')

  setting_args = group.add_argument_group('Settings')
  setting_args.add_argument("--split", dest="split", help="Split output per timestamp", action='store_true')

def fastqc_subparser(subparsers):
  group = subparsers.add_parser('fastqc', description='Summarise FastQC output', help='Summarise FastQC output')

  input_args = group.add_argument_group('Input arguments')
  input_args.add_argument("-i", "--input-dir", dest="inputdir", help="Location of assembly working directory and glob all summaries", metavar="HEAD DIR", default='~/assembly', type=str)
  input_args.add_argument("--selected-summaries", dest="selected_summaries", help="List selected summaries to analyse", metavar="SUMMARIES", type=str, nargs='+')

  setting_args = group.add_argument_group('Settings')
  setting_args.add_argument("--split", dest="split", help="Split output per timestamp", action='store_true')

def plot_subparser(subparsers):
  group = subparsers.add_parser('plot', description='Plot basic statistics', help='Plot basic statistics')

  input_args = group.add_argument_group('Input arguments')
  input_args.add_argument("-i", "--input-dir", dest="inputdir", help="Location of assembly working directory and glob all summaries", metavar="HEAD DIR", default='~/assembly', type=str)
  input_args.add_argument("--selected-summaries", dest="selected_summaries", help="List selected summaries to analyse", metavar="SUMMARIES", type=str, nargs='+')

  setting_args = group.add_argument_group('Settings')
  setting_args.add_argument("-c", "--columns", dest="columns", help="Columns to plot", metavar="COLUMNS", choices=['depth', 'N50', 'assembly_size', 'nr_reads_top_species', 'pct_top_species', 'ST', 'nr_contigs', 'all'], required=True, type=str, nargs='+')
  setting_args.add_argument("-s", "--species", dest='species', help='Species to plot', metavar="SPECIES", choices=['Escherichia coli', 'Neisseria meningitidis', 'Haemophilus influenzae'], required=True)
  setting_args.add_argument("--config", dest='config', help='Config file for plotting settings', metavar="CONFIG FILE", default='workflow/config/plot_settings.json')

def combinetyping_subparser(subparsers):
  group = subparsers.add_parser('combinetyping', description='Combine typing summaries into single table', help='Combine typing summaries')

  input_args = group.add_argument_group('Input arguments')
  input_args.add_argument("-i", "--input-dir", dest="inputdir", help="Location of assembly working directory and glob all summaries", metavar="HEAD DIR", default='~/assembly', type=str)

  setting_args = group.add_argument_group('Settings')
  setting_args.add_argument("-s", "--species", dest='species', help='Species to plot', metavar="SPECIES", choices=['Ecoli', 'Nmen'], required=True)

  output_args = group.add_argument_group('Output')
  setting_args.add_argument("-o", "--output", dest='output', help='Output table file', metavar="OUTPUT", required=True)


if __name__ == '__main__':
  main()
