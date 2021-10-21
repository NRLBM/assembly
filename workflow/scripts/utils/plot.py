#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from .read import glob_all_summaries, get_settings, read_summaries, select_columns

def plot(args):
  if args.selected_summaries is None:
    list_all_summaries = glob_all_summaries(args.inputdir)
    df = read_summaries(list_all_summaries)
  else:
    df = read_summaries(args.selected_summaries)
  settings_dict = get_settings(args.config)
  columns = select_columns(args.columns, settings_dict)
  for column in columns:
    if settings_dict[column]["type"] == "histplot":
      histplot_data(args.species, df, column, settings_dict)
    elif settings_dict[column]["type"] == "countplot":
      countplot_data(args.species, df, column, settings_dict)

def histplot_data(species, dataframe, column, settings_dict):
  dataframe_species = dataframe.query('`Top species` == @species')
  sns.set_style('whitegrid')
  sns.histplot(data = dataframe_species, x = settings_dict[column]["name"], stat='count', binwidth=settings_dict[column]['binwidth'])
  plt.title(' '.join([settings_dict[column]['name'], 'for', species]))
  plt.show()

def countplot_data(species, dataframe, column, settings_dict):
  dataframe_species = dataframe.query('`Top species` == @species')
  sns.set_style('whitegrid')
  column_name = settings_dict[column]["name"]
  sns.countplot(data = dataframe_species, x = column_name, order = dataframe_species[column_name].value_counts().index)
  plt.xticks(rotation=45)
  plt.title(' '.join([settings_dict[column]['name'], 'for', species]))
  plt.show()
