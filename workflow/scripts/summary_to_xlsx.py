#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='Convert .csv file to .xlsx file')

parser.add_argument('input', help='Input csv file', type=str)
parser.add_argument('output', help='Output xlsx file', type=str)
parser.add_argument('-d', help='Input file delimiter', type=str, default=',')

args = parser.parse_args()

import openpyxl
import pandas as pd

# Read input csv table
df = pd.read_csv(args.input, sep=args.d)

# Output to xlsx using openpyxl
df.to_excel(args.output, engine='openpyxl', index=None)
