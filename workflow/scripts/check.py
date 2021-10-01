#!/usr/bin/env python3

# Load libraries
import os
import re
import subprocess

def list_read_files():
  '''
  Find all files in input folder

  Uses a subprocess to find all R1 and R2 files as they come from the MiSeq sequencer.

  Returns
  -------
  list_files : list
    List of R1 and R2 files identified through find
  
  '''
  completed_find = subprocess.run(['find input/ -size +0 -type f -regextype sed -regex ".*/*_L001_R[12]_001.fastq.gz" -exec basename {} \;'], capture_output=True, shell=True)
  list_files = completed_find.stdout.decode("utf-8").rstrip('\n').split('\n')
  return list_files

def check_two_files(list_files):
  '''
  Check whether all samples have two files attached to them

  For all files in list_files, check whether a single R1 and a single R2 file is present.

  Parameters
  ----------
  list_files : list
    List of R1 and R2 files identified through find

  Returns
  -------
  count_dict : dictionary
    Dictionary containing sample names as keys and number of files (either R1 or R2) as integer value

  Raises
  ------
  ValueError
    If any key in the dictionary has a value which is not 2

  '''
  count_dict = {}

  for file_name in list_files:
    sample_name = re.sub('_L001_R[12]_001.fastq.gz', '', file_name)
    if (sample_name in count_dict):
      count_dict[sample_name] += 1
    else:
      count_dict[sample_name] = 1

  for key, value in count_dict.items():
    if value != 2:
      errstring = ''.join(['Sample ', str(key), ' has ', str(value), ' files in input folder where two are expected. Exiting now.'])
      raise ValueError(errstring)

  return count_dict

def check_R1_R2_per_sample(list_files, count_dict):
  '''
  Check whether every sample has an R1 and R2 file associated with them

  Goes through the count_dict dictionary and checks whether every sample is associated to exactly one R1 and one R2 file.

  Parameters
  ----------
  list_files : list
    List of R1 and R2 files identified through find
  count_dict : dictionary
    Dictionary containing sample names as keys and number of files (either R1 or R2) as integer value

  Returns
  -------
  Nothing

  Raises
  ------
  ValueError
    If a sample has not exactly one R1 and one R2 file associated with them

  '''
  count_dict_copy = count_dict.copy()

  for key, value in count_dict.items():
    R1_string = ''.join([key, '_L001_R1_001.fastq.gz'])
    R2_string = ''.join([key, '_L001_R2_001.fastq.gz'])
    if R1_string in list_files:
      count_dict_copy[key] -= 1
    if R2_string in list_files:
      count_dict_copy[key] -= 1

  for key, value in count_dict_copy.items():
    if value != 0:
      errstring = ''.join(['Sample ', key, ' does not have exactly one R1 file and exactly one R2 file. Exiting now.'])
      raise ValueError(errstring)

def check_other_input_files(list_files):
  '''
  Check whether any other files are present in input folder

  Finds all files in input folder using the os library.

  Parameters
  ----------
  list_files : list
    List of R1 and R2 files identified through find

  Returns
  -------
  Nothing

  Raises
  ------
  ValueError
    If a file is found in the input folder but not in list_files

  '''
  all_files_inputdir = os.listdir('input/')
  for file in all_files_inputdir:
    if file not in list_files:
      errstring = ''.join([file, ' is present in input dir and does not seem to be a non-empty sequencing read file. Exiting now.'])
      raise ValueError(errstring)

def calculate_pipeline_time(count_dict):
  '''
  Calculates time needed for pipeline

  Parameters
  ----------
  count_dict : dictionary
    Dictionary containing sample names as keys and number of files (either R1 or R2) as integer value

  Returns
  -------
  nr_minutes : int
    Number of minutes needed for pipeline to execute

  '''
  nr_minutes = len(count_dict) * 20
  return nr_minutes

def main():
  # List all R1 and R2 files in input folder
  list_files = list_read_files()
  # Check whether two files per sample are present
  count_dict = check_two_files(list_files)
  # Check whether R1 and R2 are present per sample
  check_R1_R2_per_sample(list_files, count_dict)
  # Check whether other input files are present
  check_other_input_files(list_files)
  # Calculate number of minutes estimated for pipeline execution and print this
  nr_minutes = calculate_pipeline_time(count_dict)
  print(nr_minutes)

if __name__ == '__main__':
  main()
