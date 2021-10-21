#!/usr/bin/env python3

import pandas as pd
from collections import Counter

from .read import glob_all_summaries, read_summary, read_summaries

def fastqc(args):
  if args.selected_summaries is None:
    list_all_summaries = glob_all_summaries(args.inputdir)
  else:
    list_all_summaries = args.selected_summaries

  df = read_summaries(list_all_summaries)

  if args.split:
    for summary in list_all_summaries:
      df = read_summary(summary)
      print(summary)
      print_fastqc_output(df)
      print('\n')
  else:
    df = read_summaries(list_all_summaries)
    print_fastqc_output(df)


def print_fastqc_output(df):
  df_fastqc = pd.DataFrame()
  list_fastqc_columns = ['FastQC warnings (pre-trim)', 'FastQC failures (pre-trim)', 'FastQC warnings (post-trim)', 'FastQC failures (post-trim)']
  for column in list_fastqc_columns:
    tmp_list = []
    df[column].apply(lambda x: [tmp_list.append(x) for x in str(x).split(';')])
    tmp_dict = Counter(tmp_list)
    tmp_dict['name'] = column
    df_fastqc = df_fastqc.append(tmp_dict, ignore_index=True)
  if 'nan' in df_fastqc.index:
    df_fastqc = df_fastqc.drop('nan', axis=1)
  df_fastqc = df_fastqc.fillna(0).transpose()
  df_fastqc.columns = df_fastqc.loc['name',:]
  df_fastqc = df_fastqc.drop('name')
  print(df_fastqc.to_markdown())
