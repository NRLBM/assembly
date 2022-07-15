#!/usr/bin/env python3

import argparse

# Three things to do
## Copy Nmen genomes to folder
## Parse metadata to suitable format
## Provide list of strain names and assembly file names

parser = argparse.ArgumentParser(description='Prepare folder for PubMLST upload.')

parser.add_argument('-m', '--metadata', help="Input folder with metadata files", dest='metadata', default='metadata_pubmlst')
parser.add_argument('--qc', help="QC report file", dest='qc', required=True, type=str)
parser.add_argument('-i', '--input', help='Input folder with genome assemblies', dest='input', required=True)
parser.add_argument('-o', '--output', help="Output folder", dest='output', required=True)
parser.add_argument('-v', '--verbose', dest='verbose', action='count', help='Increase verbosity (nothing for WARNING, -v for INFO, -vv for DEBUG)')
parser.add_argument('--force-copy', dest='force_copy', action='store_true', help='Force copying of genomes file (default: do not copy if file exists and size > 0)')
parser.add_argument('--source-values', help="TSV file without header translating NRLBM source to PubMLST allowed values (default: metadata_pubmlst/source_values.tsv", dest='source_values', default='metadata_pubmlst/source_values.tsv')

args = parser.parse_args()

import logging

if not args.verbose:
  # No increased verbosity
  logging.basicConfig(level=logging.WARNING)
elif args.verbose == 1:
  # Single -v included on command line
  logging.basicConfig(level=logging.INFO)
else:
  # -vv included on command line
  logging.basicConfig(level=logging.DEBUG)

import pandas as pd
import re
import datetime
import os
from pathlib import Path
import shutil
import sys

pd.options.mode.chained_assignment = None # SettingWithCopyWarning is turned off here to limit warning messages in logfile.txt to important warnings

logging.debug(f"Reading QC report file from {args.qc}")
df_qc = pd.read_csv(args.qc, sep='\t')
logging.debug(f"Getting all isolates from QC report into list_isolates_pass")
list_isolates_pass = df_qc[df_qc['QC_decision'] == 'PASS']['Isolate'].to_list()
logging.debug(f"list_isolates_pass contains {len(list_isolates_pass)} elements")

dict_isolates = {}

def parse_max_year():
  logging.debug(f"Parsing maximum year to filter for using datetime")
  max_year_full = str(datetime.date.today().year + 1)
  logging.debug(f"Max year to filter for is {max_year_full}")
  max_year_short = max_year_full[0] + max_year_full[2:]
  logging.debug(f"max_year_short: {max_year_short}")
  return max_year_short

def parse_onsnummer_isolate(isolate, max_year):
  logging.debug(f"Starting parse_onsnummer_isolate with match = False")
  match = match_with_suffix = match_without_suffix = match_isolate = match_onsnummer = False
  pattern_with_suffix = '2[0-9]{6}-?I{1,3}'
  pattern_without_suffix = '2[0-9]{6}'
  pattern_wrong_suffix='2[0-9]{6}-[123]'
  isolate=str(isolate)
  logging.debug(f"Searching isolate {isolate} with pattern_wrong_suffix {pattern_wrong_suffix}")
  match_wrong_suffix = re.search(pattern_wrong_suffix, isolate)
  if match_wrong_suffix:
    logging.warning(f"isolate {isolate} matches a wrong suffix (-1, -2 or -3 after a valid onsnummer).")
    logging.warning(f"TAKE EXTRA CARE WHILE SUBMITTING AS THIS NAME IS PROBABLY PARSED INCORRECTLY.")
  logging.debug(f"Searching isolate {isolate} with pattern_with_suffix {pattern_with_suffix}")
  match_with_suffix = re.search(pattern_with_suffix, isolate)
  logging.debug(f"Checking whether match_with_suffix is True")
  if match_with_suffix:
    logging.debug(f"match_with_suffix is True, saving first match to variable match_isolate")
    match_isolate = match_with_suffix.group(0)
    logging.debug(f"match_isolate is now {match_isolate}")
  else:
    logging.debug(f"match_with_suffix is False, searching isolates with pattern_without_suffix {pattern_without_suffix}")
    match_without_suffix = re.search(pattern_without_suffix, isolate)
    logging.debug(f"Checking whether match_without_suffix is True")
    if match_without_suffix:
      logging.debug(f"match_without_suffix is Truem saving first match to variable match_isolate")
      match_isolate = match_without_suffix.group(0)
      logging.debug(f"match_isolate is now {match_isolate}")
  if match_isolate:
    logging.debug(f"Searching onsnummer from match_isolate")
    match_onsnummer = re.search(pattern_without_suffix, match_isolate).group(0)
  else:
    logging.warning(f"Could not identify an onsnummer for isolate {isolate}. This isolate will NOT be included in the final output.")
  if match_onsnummer:
    logging.debug(f"Checking length of match {match_onsnummer} for isolate {isolate}")
    if len(match_onsnummer) == 7:
      logging.debug(f"Checking max_year for match {match_onsnummer} for isolate {isolate}")
      if (int(match_onsnummer[:3]) >= 200) & (int(match_onsnummer[:3]) <= int(max_year)):
        logging.debug(f"onsnummer {match_onsnummer} and isolate {match_isolate} are valid, returning this now")
        return match_isolate, match_onsnummer
  else:
    logging.warning(f"Could not identify an onsnummer for isolate {isolate}. This isolate will NOT be included in the final output.")
    return None, None

def create_nested_dirs(output):
  # Creation of multiple subdirectories can be implemented by adding name of subdir to list used by for loop
  logging.debug(f"Making directories for output")
  for subdir in ['genomes']:
    logging.debug(f"Making subdir {subdir}")
    subdir_path = Path(output).joinpath(subdir)
    logging.debug(f"Using path for subdir {subdir}: {subdir_path}")
    Path(subdir_path).mkdir(parents=True, exist_ok=True)

# Year;ons_nummer;leeftijdjr;geslacht;1;2;4;microorganisme;Typering - serogroep

def filter_metadata(metadata_path, isolate_dict):
  logging.debug(f"Reading df from metadata_path {metadata_path}")
  df = pd.read_csv(metadata_path, sep=';', dtype=str, encoding='latin_1', encoding_errors='replace')
  logging.debug(f"Creating new df to write results to")
  new_df = pd.DataFrame(columns=df.columns)
  logging.debug(f"Iterating over items in isolate_dict {isolate_dict}")
  for key, value in isolate_dict.items():
    logging.debug(f"Finding entries in df where ons_nummer matches {value[1]}")
    tmp_df = df[df['ons_nummer'] == value[1]]
    logging.debug(f"Concatenating df of shape {tmp_df.shape} to new df with results")
    new_df = pd.concat([new_df, tmp_df])
  logging.debug(f"Finding duplicates and reporting these to warnings")
  duplicates = new_df[new_df.duplicated(keep=False)]
  for index, row in duplicates.iterrows():
    logging.warning(f"DUPLICATE: {row.to_dict()}")
  logging.debug(f"Dropping {duplicates.shape[0]/2} duplicates")
  new_df = new_df.drop_duplicates(ignore_index=True)
  logging.debug(f"Returning new df with shape {new_df.shape}")
  return new_df

def get_source_values(source_path):
  dict = {}
  logging.debug(f"Opening file containing source conversion values using source path {source_path}")
  with open(source_path, 'r+') as file:
    logging.debug(f"Calling readlines on file")
    lines = file.readlines()
  logging.debug(f"Looping through lines")
  for line in lines:
    logging.debug(f"Parsing key and value from line {line}")
    key, value = line.split('\t')
    stripped_value = value.rstrip('\n')
    logging.debug(f"Adding value {stripped_value} to dict for key {key}")
    dict[key] = stripped_value
  logging.debug(f"Returning dict of size {len(dict)}")
  return dict

def parse_metadata(df, source_dict):
  logging.debug(f"Setting up dicts to convert values")
  geslacht_dict = {'M': 'male', 'V': 'female'}
  microorganisme_dict = {'Neisseria meningitidis': 'Neisseria meningitidis'}
  serogroep_dict = {'A': 'A', 'B': 'B', 'C': 'C', 'E': 'E', '29E': 'E', 'W': 'W', 'W135': 'W',
                'X': 'X', 'Y': 'Y', 'Z': 'Z', 'Niet Groepeerbaar': 'NG', 'H': 'H', 'J': 'J',
                'K': 'K', 'L': 'L', 'NG': 'NG'}
  logging.debug(f"Creating a list of necessary columns")
  list_necessary_columns = ['ons_nummer', 'Year', 'geslacht', '1', '2', 'microorganisme', 'Typering - serogroep']
  logging.debug(f"Looping through columns to assert these are present")
  for col in list_necessary_columns:
    assert col in df.columns, f"Column {col} is not found in metadata table, please check the table and try again"
  logging.debug(f"Checking if source column '3' is in df.columns")
  if '3' in df.columns:
    logging.debug(f"Checking if source column '4' is in df.columns")
    if '4' in df.columns:
      logging.debug(f"Concatenating columns 1, 2, 3 and 4 for source")
      df['source_raw'] = df['1'].fillna('') + df['2'].fillna('') + df['3'].fillna('') + df['4'].fillna('')
    else:
      logging.debug(f"Concatenating columns 1, 2 and 3 for source")
      df['source_raw'] = df['1'].fillna('') + df['2'].fillna('') + df['3'].fillna('')
  else:
    logging.debug(f"Checking if source column '4' is in df.columns")
    if '4' in df.columns:
      logging.debug(f"Concatenating columns 1, 2 and 4 for source")
      df['source_raw'] = df['1'].fillna('') + df['2'].fillna('') + df['4'].fillna('')
    else:
      logging.debug(f"Concatenating columns 1 and 2 for source")
      df['source_raw'] = df['1'].fillna('') + df['2'].fillna('')
  logging.debug(f"Constructing filter for samples which only have Liquor as source")
  filter_liquor_only = df['source_raw'] == 'Liquor'
  logging.debug(f"Checking if samples have only Liquor as source")
  if sum(filter_liquor_only) > 0:
    logging.warning(f"A total of {sum(filter_liquor_only)} entries were excluded because these were listed as only liquor")
    logging.warning(f"Note that other entries for the same ons_nummer can still be present in output")
    df = df[~filter_liquor_only]
  logging.debug(f"Mapping sources for columns 1, 2, 3 and 4")
  df['source_1'] = df['1'].map(source_dict).fillna('')
  df['source_2'] = df['2'].map(source_dict).fillna('')
  if '3' in df.columns:
    df['source_3'] = df['3'].map(source_dict).fillna('')
  if '4' in df.columns:
    df['source_4'] = df['4'].map(source_dict).fillna('')
  logging.debug(f"Mapping serogroup")
  df['serogroup'] = df['Typering - serogroep'].map(serogroep_dict).fillna('')
  logging.debug(f"Mapping species")
  df['species'] = df['microorganisme'].map(microorganisme_dict).fillna('Neisseria sp.')
  logging.debug(f"Mapping sex")
  df['sex'] = df['geslacht'].map(geslacht_dict).fillna('')
  logging.debug(f"Adding The Netherlands as country")
  df['country'] = 'The Netherlands'
  logging.debug(f"Renaming ons_nummer to isolate")
  df['isolate'] = df['ons_nummer']
  logging.debug(f"Renaming Year to year")
  df['year'] = df['Year']
  logging.debug(f"Adding Illumina as method")
  df['sequence_method'] = 'Illumina'
  logging.debug(f"Returning df with {df.shape}")
  return df

def parse_source_suffix(row):
  logging.debug(f"Searching III in isolate name")
  isolate = str(row['isolate'])
  match_other = re.search('III', isolate)
  if match_other:
    logging.warning(f"Isolate {isolate} was not isolated from blood/CSF.")
    logging.warning(f"Please translate the source and add this to the comments field when uploading.")
    logging.debug(f"III was found in isolate name, returning source_3 column value")
    return row['source_3']
  else:
    logging.debug(f"Searching II in isolate name")
    match_blood = re.search('II', row['isolate'])
    if match_blood:
      logging.debug(f"II was found in isolate name, returning source_2 column value")
      return row['source_2']
    else:
      logging.debug(f"Searching for I in isolate name, returning source_1 column value")
      match_CSF = re.search('I', row['isolate'])
      if match_CSF:
        logging.debug(f"I was found in isolate name, returning source_1 column value")
        return row['source_1']

def finalise_metadata(df, dict_isolates):
  logging.debug(f"Defining combine columns")
  combine_cols = ['source_1', 'source_2']
  skip_multiple_source = skip_single_source = skip_no_source = False
  df = df.copy(deep=True)

  logging.debug(f"Checking if source_3 is present")
  if 'source_3' in df.columns:
    combine_cols.append('source_3')
  logging.debug(f"Checking if source_4 is present")
  if 'source_4' in df.columns:
    combine_cols.append('source_4')
  logging.debug(f"Check how many sources have been listed for this entry using combine_cols: {combine_cols}")
  df['sum_sources'] = (df[combine_cols] != '').sum(axis=1)

  # Multiple sources
  logging.debug(f"Finding entries which have more than one source, saving these in new df")
  df_multiple_source = df[df['sum_sources'] > 1]
  if len(df_multiple_source) > 0:
    logging.debug(f"Getting onsnummers from new df with multiple sources")
    onsnummers = list(df_multiple_source['ons_nummer'])
    logging.debug(f"Creating list_isolates_mult to store isolate names and list_onsnummers_mult for onsnummers")
    list_isolates_mult = []
    list_onsnummers_mult = []
    logging.debug(f"Iterating over dict_isolates")
    for key, value in dict_isolates.items():
      logging.debug(f"Checking if value[1]: {value[1]} is in onsnummers")
      if value[1] in onsnummers:
        logging.debug(f"value[1] was found in list_onsnummers")
        logging.debug(f"value[0]: {value[0]} is appended to list isolates")
        logging.debug(f"value[1]: {value[1]} is appended to list onsnummers")
        list_isolates_mult.append(value[0])
        list_onsnummers_mult.append(value[1])
    logging.debug(f"Dropping isolate column from df_multiple_source to prevent collision in later merge")
    df_multiple_source = df_multiple_source.drop('isolate', axis=1)
    logging.debug(f"Create new_df_multiple_source from list_isolates_mult and list_onsnummers_mult")
    new_df_multiple_source = pd.DataFrame.from_dict({'isolate': list_isolates_mult, 'ons_nummer': list_onsnummers_mult}, dtype='str')
    logging.debug(f"Merge new_df_multiple_source with df_multiple_source on onsnummer")
    new_df_multiple_source = new_df_multiple_source.merge(df_multiple_source, on='ons_nummer', how='left')
    logging.debug(f"Parse source from isolate names using function parse_source_suffix")
    new_df_multiple_source.loc[:,'source'] = new_df_multiple_source.apply(parse_source_suffix, axis=1)
  else:
    logging.debug(f"No onsnummers were listed as from multiple sources")
    skip_multiple_source = True

  # Single source
  logging.debug(f"Saving entries with single source to df_single_source")
  df_single_source = df[df['sum_sources'] == 1]
  if len(df_single_source) > 0:
    logging.debug(f"Concatenating combine_cols {combine_cols} into source column")
    df_single_source['source'] = df_single_source[combine_cols].apply(lambda x: ''.join(x.values.astype(str)), axis=1)
  else:
    logging.debug(f"No onsnummers were listed as from a single source")
    skip_single_source = True

  # No source
  logging.debug(f"Saving entries without source to df_no_source")
  df_no_source = df[df['sum_sources'] == 0]
  if len(df_no_source) > 0:
    logging.debug(f"Adding missing empty string to source column")
    df_no_source.loc[:,'source'] = ''
  else:
    logging.debug(f"No onsnummers were listed as without a source")
    skip_no_source = True

  if all([skip_multiple_source, skip_single_source, skip_no_source]):
    logging.warning(f"No samples were found. Aborting...")
    sys.exit(1)
  elif all([skip_multiple_source, skip_single_source]):
    total_df = df_no_source
  elif all([skip_single_source, skip_no_source]):
    total_df = new_df_multiple_source
  elif all([skip_multiple_source, skip_no_source]):
    total_df = df_single_source
  elif skip_multiple_source:
    total_df = pd.concat([df_single_source, df_no_source])
  elif skip_single_source:
    total_df = pd.concat([new_df_multiple_source, df_no_source])
  elif skip_no_source:
    total_df = pd.concat([new_df_multiple_source, df_single_source])
  else:
    logging.debug(f"Concatenating all three dfs into total_df")
    total_df = pd.concat([new_df_multiple_source, df_single_source, df_no_source])

  logging.debug(f"Parse assembly_filename from isolate column")
  total_df['assembly_filename'] = total_df['isolate'] + '.fasta'
  logging.debug(f"Parsing disease from source (necessary for PubMLST)")
  total_df['disease'] = total_df['source'].apply(lambda x: 'invasive (unspecified/other)' if x in ['blood', 'CSF', 'joint fluid'] else '')

  return total_df

def create_upload_sheet(df, isolate_dict):
  logging.debug(f"Create tuple of columns needed for upload sheet")
  needed_columns = ('isolate', 'aliases', 'references', 'country', 'region', 'town_or_city', 'year', 'date_sampled',
		'date_received', 'non_culture', 'epidemiological_year', 'age_yr', 'age_range', 'age_mth', 'sex', 'disease',
		'source', 'epidemiology', 'species', 'serogroup', 'genogroup', 'serotype', 'sero_subtype', 'amoxicillin_mic_sign',
		'amoxicillin_mic', 'azithromycin_mic_sign', 'azithromycin_mic', 'cefixime_mic_sign', 'cefixime_mic',
		'cefotaxime_mic_sign', 'cefotaxime_mic', 'ceftriaxone_mic_sign', 'ceftriaxone_mic', 'chloramphenicol_mic_sign',
		'chloramphenicol_mic', 'ciprofloxacin_mic_sign', 'ciprofloxacin_mic', 'penicillin_mic_sign', 'penicillin_mic',
		'rifampicin_mic_sign', 'rifampicin_mic', 'spectinomycin_mic_sign', 'spectinomycin_mic', 'sulphonamide_mic_sign',
		'sulphonamide_mic', 'tetracycline_mic_sign', 'tetracycline_mic', 'bioproject_accession', 'biosample_accession',
		'ENA_run_accession', 'comments', 'assembly_filename', 'sequence_method')
  logging.debug(f"Looping through needed_columns: {needed_columns}")
  for col in needed_columns:
    logging.debug(f"Checking whether column {col} is already present in df")
    if col not in df.columns:
      logging.debug(f"column {col} is added as an empty column")
      df[col] = ''
    else:
      logging.debug(f"columns {col} was already present")
  logging.debug(f"Creating cleaned_df from df containing all data")
  logging.debug(f"Shape of df is {df.shape}")
  cleaned_df = df[list(needed_columns)]
  logging.debug(f"Shape of cleaned_df is {cleaned_df.shape}")
  return cleaned_df

def save_metadata_sheet(df, output_dir):
  logging.debug(f"Creating output path from output_dir: {output_dir} and PubMLST_upload_sheet.txt")
  output_path = Path(output_dir).joinpath('PubMLST_upload_sheet.txt')
  logging.debug(f"Saving output sheet to output_path: {output_path}")
  df.to_csv(output_path, sep='\t', index=False)

def copy_genomes(dict_isolates, input, output, force_copy):
  logging.debug(f"Parsing genome_dir from {output} and 'genomes'")
  genome_dir = Path(output).joinpath('genomes')
  logging.debug(f"Asserting genome_dir {genome_dir} exists")
  assert os.path.exists(genome_dir), f"genome_dir {genome_dir} does not seem to exist. Was function create_nested_dirs successful?"
  for raw_isolate, clean_isolate in dict_isolates.items():
    copy_flag = False
    logging.debug(f"Constructing clean_isolate_path for isolate {isolate}")
    clean_isolate_path = Path(genome_dir).joinpath(clean_isolate[0] + '.fasta')
    logging.debug(f"Checking if clean_isolate_path {clean_isolate_path} already exists and is not an empty file")
    clean_isolate_exists = os.path.exists(clean_isolate_path)
    logging.debug(f"If clean_isolate_path: {clean_isolate_path} exists, check size")
    if clean_isolate_exists:
      clean_isolate_size = os.path.getsize(clean_isolate_path)
    logging.debug(f"If clean_isolate_path: {clean_isolate_path} does not exist, copy genome file")
    if not clean_isolate_exists:
      copy_flag = True
    else:
      if clean_isolate_size == 0:
        logging.debug(f"If clean_isolate_path: {clean_isolate_path} has size 0, copy genome file")
        copy_flag = True
    if force_copy == True:
      logging.debug(f"--force-copy argument is used, copy genome file")
      copy_flag = True
    logging.debug(f"Check if copy_flag: {copy_flag} is True")
    if copy_flag:
      logging.debug(f"Constructing raw_isolate_path for isolate {isolate}")
      raw_isolate_path = Path(input).joinpath(raw_isolate + '.fasta')
      logging.debug(f"Copying isolate {isolate} from {raw_isolate_path} to {clean_isolate_path}")
      shutil.copy2(raw_isolate_path, clean_isolate_path)

## Running functions
logging.debug(f"Calling function parse_max_year")
max_year_short = parse_max_year()

logging.debug(f"Looping through list_isolates_pass: {list_isolates_pass}")
for isolate in list_isolates_pass:
  logging.debug(f"Calling function parse_onsnummer_isolate on isolate {isolate}")
  isolate_parsed, onsnummer_parsed = parse_onsnummer_isolate(isolate, max_year_short)
  if (isolate_parsed is not None) & (onsnummer_parsed is not None):
    dict_isolates[isolate] = (isolate_parsed, onsnummer_parsed)

logging.debug(f"Calling function create_nested_dirs on output folder {args.output}")
create_nested_dirs(args.output)

source_dict = get_source_values(args.source_values)

df = filter_metadata(args.metadata, dict_isolates)
df_parsed = parse_metadata(df, source_dict)
df_parsed2 = finalise_metadata(df_parsed, dict_isolates)
complete_df = create_upload_sheet(df_parsed2, dict_isolates)
save_metadata_sheet(complete_df, args.output)
copy_genomes(dict_isolates, args.input, args.output, args.force_copy)

# Write all doubts/warnings to warnings.txt
