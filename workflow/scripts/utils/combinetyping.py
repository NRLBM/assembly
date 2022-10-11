#!/usr/bin/env python3

from .parseyear import parse_year
from .read import glob_all_typing_summaries, read_typing_summaries
import pandas as pd

def combinetyping(args):
  if args.species == "Escherichia coli":
    args.species = "Ecoli"
  elif args.species == "Neisseria meningtidis":
    args.species == "Nmen"
  elif args.species == "Streptococcus pyogenes":
    args.species == "Spyo"
  list_all_summaries = glob_all_typing_summaries(args.inputdir, args.species)
  df = read_typing_summaries(list_all_summaries)
  year_series = df['Isolate'].apply(parse_year)
  df.insert(1, 'year', year_series)
  df.to_csv(args.output, sep = '\t', index=False)
