#!/usr/bin/env python3

__author__ = "Boas van der Putten"
__version__ = "0.0.1"

import argparse

parser = argparse.ArgumentParser(description='Convert .csv file to .xlsx file')

parser.add_argument('input', help='Input csv file', type=str)
parser.add_argument('versions', help='Versions file', type=str)
parser.add_argument('output', help='Output xlsx file', type=str)
parser.add_argument('-d', help='Input file delimiter', type=str, default=',')
parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))

args = parser.parse_args()

import openpyxl
import pandas as pd

# Read input csv table
df = pd.read_csv(args.input, sep=args.d)
df_versions = pd.read_csv(args.versions, sep=",")

# Output to xlsx using openpyxl
with pd.ExcelWriter(args.output, engine='openpyxl') as writer:
  df.to_excel(writer, index=None, sheet_name='data')
  df_versions.to_excel(writer, index=None, sheet_name='versions')
