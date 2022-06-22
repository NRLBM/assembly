#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='Find and copy input file to output file, based on expected pattern.')

parser.add_argument('--input', dest='input', required=True, type=str, help='Input directory')
parser.add_argument('--output', dest='output', required=True, type=str, nargs='+', help='Output files')
parser.add_argument('--pattern', dest='pattern', type=str, default='microbesng', choices=['microbesng'], help='Expected pattern of input file. Default: MicrobesNG. More patterns will be added in the future.')

parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="Increase verbosity (nothing for WARNING, -v for INFO, -vv for DEBUG)")

args = parser.parse_args()

import os
import sys
import subprocess
import logging

# Script uses f-strings, check whether Python>=3.6
assert sys.version_info >= (3, 7), "This script requires at least Python version 3.7"

if not args.verbose:
  # No increased verbosity
  logging.basicConfig(level=logging.WARNING)
elif args.verbose == 1:
  # Single -v included on command line
  logging.basicConfig(level=logging.INFO)
else:
  # -vv included on command line
  logging.basicConfig(level=logging.DEBUG)

# Decide settings based on pattern
if args.pattern == 'microbesng':
  logging.debug(f"Deciding settings")
  regex_pattern = '[0-9]\{6\}_'
  suffix = '_trimmed.fastq.gz'
  logging.debug(f"Settings for {args.pattern}:")
  logging.debug(f"regex_pattern: {regex_pattern}")
  logging.debug(f"suffix: {suffix}")

def check_output(list_output, pattern):
  '''
  Check whether two valid output file names are provided.

  Parameters
  ----------
  list_output : list
    List of two valid output files
  pattern: str
    Pattern settings to use

  Raises
  ------
  AssertionError
    If list of output files does not contain exactly two elements
  AssertionError
    If first output file does not end with "_1.fastq.gz"
  AssertionError
    If second output file does not end with "_2.fastq.gz"
  AssertionError
    If basenames of fw and rv output files do not match

  '''
  logging.debug(f"Asserting two output files are provided")
  assert len(list_output) == 2, "Output file list is not length 2"
  logging.debug(f"Unpacking output files into output_fw and output_rv")
  output_fw, output_rv = list_output
  if pattern == 'microbesng':
    logging.debug(f"Asserting output_fw ({output_fw}) ends with \"_1.fastq.gz\"")
    assert output_fw[-11:] == '_1.fastq.gz', f"First file of {list_output} does not end in \"_1.fastq.gz\""
    logging.debug(f"Asserting output_rv ({output_rv}) ends with \"_2.fastq.gz\"")
    assert output_rv[-11:] == '_2.fastq.gz', f"Second file of {list_output} does not end in \"_2.fastq.gz\""
    output_fw_base = os.path.splitext(os.path.splitext(os.path.basename(output_fw))[0])[0]
    output_rv_base = os.path.splitext(os.path.splitext(os.path.basename(output_rv))[0])[0]
    assert output_fw_base[:-1] == output_rv_base[:-1], f"Basenames of fw and rv output files does not match: {output_fw_base} and {output_rv_base}"

def find_inputfile(output_file, regex_pattern, suffix, inputdir):
  '''
  Find input file using unix find.

  Parameters
  ----------
  output_file : str
    Path to output file to find input file for
  regex_pattern : str
    Regex pattern to include in find command
  suffix : str
    Expected suffix of input file
  inputdir : str
    Path to directory where input files must be searched

  Returns
  -------
  input_file : str
    Path to input file identified by unix find

  Raises
  ------
  AssertionError
    If no files are found using unix find
  AssertionError
    If multiple files are found using unix find

  '''
  logging.debug(f"Finding basename for output_file {output_file}")
  output_file_base = os.path.splitext(os.path.splitext(os.path.basename(output_file))[0])[0]
  logging.debug(f"Constructing input pattern to find based on inputdir, regex_pattern, output_file_base and suffix")
  input_pattern = ''.join([inputdir, '/', regex_pattern, output_file_base, suffix])
  logging.debug(f"Constructing find command")
  # Find does not use size >0 bytes because a gzipped empty file is more than 0 bytes. Anything up to 100 bytes is now excluded.
  find_command = ''.join(['find ', inputdir, ' -size +100c -type f -regextype sed -regex "', input_pattern, '" -print'])
  logging.info(f"Running find command: {find_command}")
  completed_find = subprocess.run(find_command, capture_output=True, shell=True)
  logging.debug(f"Decoding captured output from find_command")
  input_file = completed_find.stdout.decode("utf-8").rstrip('\n')
  assert input_file != "", f"Found no files matching find_command: {find_command}"
  assert len(input_file.split('\n')) == 1, f"Found multiple files matching find_command: {find_command}"
  logging.info(f"Found {input_file}")
  return input_file

def copy_inputfile(input_file, output_file):
  '''
  Copy input file to output file using unix cp

  Parameters
  ----------
  input_file : str
    Path to input file

  output_file : str
    Path to output file

  '''
  logging.debug(f"Constructing copy command for {input_file} and {output_file}")
  copy_command = ''.join(['cp ', input_file, ' ', output_file])
  logging.info(f"Running copy command: {copy_command}")
  subprocess.run(copy_command, shell=True)

logging.debug(f"Calling check_output function")
check_output(args.output, args.pattern)

for output_file in args.output:
  logging.debug(f"Starting {output_file} out of {args.output}")
  logging.debug(f"Calling find_inputfile function")
  input_file = find_inputfile(output_file, regex_pattern, suffix, args.input)
  logging.debug(f"Calling copy_inputfile function")
  copy_inputfile(input_file, output_file)
