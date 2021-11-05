#!/usr/bin/env python3

import pandas as pd
import json
import os
import glob

def get_settings(json_filepath):
  with open(json_filepath, 'r') as f:
    plot_settings = json.load(f)
  return plot_settings

def glob_all_summaries(inputdir):
  expanded_inputdir = os.path.expanduser(inputdir)
  glob_path = os.path.join(expanded_inputdir, 'backup/*/summary.csv')
  list_all_summaries = glob.glob(glob_path)
  return list_all_summaries

def glob_all_typing_summaries(inputdir, species):
  expanded_inputdir = os.path.expanduser(inputdir)
  glob_path = os.path.join(expanded_inputdir, ''.join(['backup/*/', species, '_typing_summary.tsv']))
  list_all_summaries = glob.glob(glob_path)
  return list_all_summaries

def read_summary(summary_file):
  df = pd.read_csv(summary_file, sep = ',', dtype = {'ST': str, 'N50': 'Int64'})
  return df

def read_summaries(list_summaries):
  df = pd.DataFrame()
  for summary_file in list_summaries:
    # print('Shape of ' + summary_file)
    tmp_df = pd.read_csv(summary_file, sep = ',', dtype = {'ST': str, 'N50': 'Int64'})
    # print(tmp_df.shape)
    df = pd.concat([df, tmp_df])
  df = df.reset_index()
  return df

def read_typing_summaries(list_summaries):
  df = pd.DataFrame()
  for summary_file in list_summaries:
    # print('Shape of ' + summary_file)
    tmp_df = pd.read_csv(summary_file, sep = '\t', dtype = str)
    # print(tmp_df.shape)
    df = pd.concat([df, tmp_df])
  df = df.reset_index().drop('index', axis=1)
  return df

def select_columns(columns, settings_dict):
  if 'all' not in columns:
    return columns
  else:
    columns = list(settings_dict.keys())
    return columns
