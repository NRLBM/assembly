#!/usr/bin/env python3

__author__ = "Boas van der Putten"
__version__ = "0.0.1"

import argparse
import zipfile

parser = argparse.ArgumentParser(description='Summarise assembled isolates')

parser.add_argument('samples', nargs='+', help="List of sample names to summarise", type=str)
parser.add_argument('--mlst', dest="mlst", help="mlst output directory", type=str, default='mlst_out')
parser.add_argument('--quast', dest="quast", help="Quast directory", type=str, default='quast_out')
parser.add_argument('--coverage', dest="coverage", help="Coverage directory", type=str, default='coverage_out')
parser.add_argument('--kraken', dest="kraken", help="Kraken2 output directory", type=str, default='kraken_out')
parser.add_argument('--fastqc-pre', dest="fastqc_pre", help="FastQC pre-trimming output directory", type=str, default='fastqc_pre_out')
parser.add_argument('--fastqc-post', dest="fastqc_post", help="FastQC post-trimming output directory", type=str, default='fastqc_post_out')
parser.add_argument('--timestamp', dest="timestamp", help="Timestamp of analysis", type=str, required=True)
parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))

args = parser.parse_args()

# Store samples as isolate_list
isolate_list = list(args.samples)

# Make output prefix based on timestamp
output_prefix = 'tmp_data/' + str(args.timestamp) + '/'

# Print header to STDOUT
print('Isolate', 'Top species', 'Pct top species', 'Number of reads top species', 'ST scheme (PubMLST)', 'ST', 'Number of contigs', 'N50', 'Largest contig', 'Total assembly size', 'Sequencing depth', 'FastQC warnings (pre-trim)', 'FastQC failures (pre-trim)', 'FastQC warnings (post-trim)', 'FastQC failures (post-trim)', sep = ',')

def customreadlines(file):
  '''
  Read all lines of a file and close connection.

  Parameters
  ----------
  file : str
    File name to open

  Returns
  -------
  tmp.readlines() : list
    List containing all lines

  '''
  tmp = open(file)
  return tmp.readlines()
  tmp.close()

# Loop through all isolates
for isolate in isolate_list:
  # Remove trailing newline
  isolate = isolate.rstrip('\n')

  # Get Kraken2 report and read lines
  lines = customreadlines(output_prefix + args.kraken + '/' + isolate + '.txt')

  # Loop through Kraken2 report
  for line in lines:
    # If taxonomic level equals species ("S")
    if line.split('\t')[3] == 'S':
      # Parse data: species, percentage or reads and absolute number of reads. Then break
      SPECIES = line.split('\t')[5].lstrip(' ').rstrip('\n')
      PCT = line.split('\t')[0].lstrip(' ')
      READS = line.split('\t')[2]
      break

  # Get ST
  t = open(output_prefix + args.mlst + '/' + isolate + ".txt")
  # Read only a single line
  line = t.readline()
  # Parse identified PubMLST scheme and which ST was identified
  ST_SCHEME = str(line.split('\t')[1])
  ST = str(line.split('\t')[2].rstrip('\n'))
  t.close()

  # Get quast report.tsv
  lines = customreadlines(output_prefix + args.quast + '/' + isolate + "/report.tsv")

  # Go through lines and check for all metrics per line. Parse if the correct metric is identified
  for line in lines:
    if line.split('\t')[0] == '# contigs':
      CONTIGS = line.split('\t')[1].rstrip('\n')
    if line.split('\t')[0] == 'Largest contig':
      LARGEST = line.split('\t')[1].rstrip('\n')
    if line.split('\t')[0] == 'Total length':
      SIZE = line.split('\t')[1].rstrip('\n')
    if line.split('\t')[0] == 'N50':
      N50 = line.split('\t')[1].rstrip('\n')

  # Get coverage, single line, single column
  c = open(output_prefix + args.coverage + '/' + isolate + '.txt')
  COVERAGE = str(c.readline().rstrip('\n'))
  c.close()

  # Get FastQC pre-trimming for both R1 and R2
  WARN_PRE = list()
  FAIL_PRE = list()

  # Do this for R1 and R2
  for read in ['1', '2']:
    # Open temporary connection to zipfile
    with zipfile.ZipFile(output_prefix + args.fastqc_pre + '/' + isolate + '_R' + read + '_fastqc.zip') as myzip:
      # Open unzipped summary file
      with myzip.open(isolate + '_' + read + '_fastqc/summary.txt') as myfile:
        # Read whole summary at once
        lines = myfile.readlines()

    # Loop through lines and decode line before checking warnings and errors
    for line_bytes in lines:
      line = line_bytes.decode("utf-8")
      if line.split('\t')[0] == 'WARN':
        # If warning is found, append to list
        WARN_PRE.append(line.split('\t')[1])
      if line.split('\t')[0] == 'FAIL':
        FAIL_PRE.append(line.split('\t')[1])

  # Join all warnings for R1 and R2 together and leave only unique ones
  WARN_PRE_FINAL = ';'.join(sorted(set(WARN_PRE)))
  # Idem for failures
  FAIL_PRE_FINAL = ';'.join(sorted(set(FAIL_PRE)))

  # Set to "NA" if no warnings or failures
  if WARN_PRE_FINAL == '':
    WARN_PRE_FINAL = 'NA'

  if FAIL_PRE_FINAL == '':
    FAIL_PRE_FINAL = 'NA'

  # Do the same for FastQC reports after trimming
  WARN_POST = list()
  FAIL_POST = list()

  for read in ['1', '2']:
    with zipfile.ZipFile(output_prefix + args.fastqc_post + '/' + isolate + '_R' + read + '_fastqc.zip') as myzip:
      with myzip.open(isolate + '_' + read + '_corrected_fastqc/summary.txt') as myfile:
        lines = myfile.readlines()

    for line_bytes in lines:
      line = line_bytes.decode("utf-8")
      if line.split('\t')[0] == 'WARN':
        WARN_POST.append(line.split('\t')[1])
      if line.split('\t')[0] == 'FAIL':
        FAIL_POST.append(line.split('\t')[1])

  WARN_POST_FINAL = ';'.join(sorted(set(WARN_POST)))
  FAIL_POST_FINAL = ';'.join(sorted(set(FAIL_POST)))

  if WARN_POST_FINAL == '':
    WARN_POST_FINAL = 'NA'

  if FAIL_POST_FINAL == '':
    FAIL_POST_FINAL = 'NA'

  # Finally print all collected variables in one line to STDOUT
  print(isolate, SPECIES, PCT, READS, ST_SCHEME, ST, CONTIGS, N50, LARGEST, SIZE, COVERAGE, WARN_PRE_FINAL, FAIL_PRE_FINAL, WARN_POST_FINAL, FAIL_POST_FINAL, sep = ',')
