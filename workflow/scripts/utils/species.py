#!/usr/bin/env python3

import pandas as pd

from .read import glob_all_summaries, read_summary, read_summaries

def species(args):
  if args.selected_summaries is None:
    list_all_summaries = glob_all_summaries(args.inputdir)
  else:
    list_all_summaries = args.selected_summaries

  df = read_summaries(list_all_summaries)

  if args.split:
    for summary in list_all_summaries:
      df = read_summary(summary)
      print(summary)
      print_species_count(df)
      print('\n')
  else:
    df = read_summaries(list_all_summaries)
    print_species_count(df)


def print_species_count(df):
  species_counts = df['Top species'].value_counts()
  df_species_counts = pd.DataFrame({'Species': species_counts.index, 'Count': species_counts})
  total_count = df_species_counts['Count'].sum()
  df_species_counts = df_species_counts.append({'Species': 'Total', 'Count': total_count}, ignore_index=True)
  print(df_species_counts.to_markdown(index=None))

