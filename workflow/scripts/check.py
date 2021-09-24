#!/usr/bin/env python3

# Load libraries
import os
import re
import subprocess

def list_read_files():
  # Find all files using a find command
  completed_find = subprocess.run(['find input/ -size +0 -type f -regextype sed -regex ".*/*_L001_R[12]_001.fastq.gz" -exec basename {} \;'], capture_output=True, shell=True)
  list_files = completed_find.stdout.decode("utf-8").rstrip('\n').split('\n')
  return list_files

def check_two_files(list_files):
  # Check whether all samples have two files attached to them
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
  # Check whether every sample which passed the previous test has a R1 and R2 associated with them
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
  # Check whether any other files are present in input folder
  all_files_inputdir = os.listdir('input/')
  for file in all_files_inputdir:
    if file not in list_files:
      errstring = ''.join([file, ' is present in input dir and does not seem to be a non-empty sequencing read file. Exiting now.'])
      raise ValueError(errstring)

def calculate_pipeline_time(count_dict):
  nr_minutes = len(count_dict) * 20
  return nr_minutes

def main():
  list_files = list_read_files()
  count_dict = check_two_files(list_files)
  check_R1_R2_per_sample(list_files, count_dict)
  check_other_input_files(list_files)
  nr_minutes = calculate_pipeline_time(count_dict)
  print(nr_minutes)

if __name__ == '__main__':
  main()
