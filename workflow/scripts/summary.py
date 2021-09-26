#!/usr/bin/env python3

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

args = parser.parse_args()

isolate_list = list(args.samples)

output_prefix = 'tmp_data/' + str(args.timestamp) + '/'

print('Sample', 'Top species', 'Pct top species', 'Number of reads top species', 'ST scheme (PubMLST)', 'ST', 'Number of contigs', 'N50', 'Largest contig', 'Total assembly size', 'Sequencing depth', 'FastQC warnings (pre-trim)', 'FastQC failures (pre-trim)', 'FastQC warnings (post-trim)', 'FastQC failures (post-trim)', sep = ',')

def customreadlines(file):
  tmp = open(file)
  return tmp.readlines()
  tmp.close()

for isolate in isolate_list:
  # Remove trailing newline
  isolate = isolate.rstrip('\n')

  # Get Kraken2
  lines = customreadlines(output_prefix + args.kraken + '/' + isolate + '.txt')

  for line in lines:
    if line.split('\t')[3] == 'S':
      SPECIES = line.split('\t')[5].lstrip(' ').rstrip('\n')
      PCT = line.split('\t')[0].lstrip(' ')
      READS = line.split('\t')[2]
      break

  # Get ST
  t = open(output_prefix + args.mlst + '/' + isolate + ".txt")
  line = t.readline()
  ST_SCHEME = str(line.split('\t')[1])
  ST = str(line.split('\t')[2])
  t.close()

  # Get quast
  lines = customreadlines(output_prefix + args.quast + '/' + isolate + "/report.tsv")

  for line in lines:
    if line.split('\t')[0] == '# contigs':
      CONTIGS = line.split('\t')[1].rstrip('\n')
    if line.split('\t')[0] == 'Largest contig':
      LARGEST = line.split('\t')[1].rstrip('\n')
    if line.split('\t')[0] == 'Total length':
      SIZE = line.split('\t')[1].rstrip('\n')
    if line.split('\t')[0] == 'N50':
      N50 = line.split('\t')[1].rstrip('\n')

  # Get coverage
  c = open(output_prefix + args.coverage + '/' + isolate + '.txt')
  COVERAGE = str(c.readline().rstrip('\n'))
  c.close()

  # Get FastQC pre-trimming for both R1 and R2
  WARN_PRE = list()
  FAIL_PRE = list()

  for read in ['R1', 'R2']:
    with zipfile.ZipFile(output_prefix + args.fastqc_pre + '/' + isolate + '_' + read + '_fastqc.zip') as myzip:
      with myzip.open(isolate + '_L001_' + read + '_001_fastqc/summary.txt') as myfile:
        lines = myfile.readlines()

    for line_bytes in lines:
      line = line_bytes.decode("utf-8")
      if line.split('\t')[0] == 'WARN':
        WARN_PRE.append(line.split('\t')[1])
      if line.split('\t')[0] == 'FAIL':
        FAIL_PRE.append(line.split('\t')[1])

  WARN_PRE_FINAL = ';'.join(sorted(set(WARN_PRE)))
  FAIL_PRE_FINAL = ';'.join(sorted(set(FAIL_PRE)))

  if WARN_PRE_FINAL == '':
    WARN_PRE_FINAL = 'NA'

  if FAIL_PRE_FINAL == '':
    FAIL_PRE_FINAL = 'NA'

  # Get FastQC post-trimming for both R1 and R2
  WARN_POST = list()
  FAIL_POST = list()

  for read in ['R1', 'R2']:
    with zipfile.ZipFile(output_prefix + args.fastqc_post + '/' + isolate + '_' + read + '_fastqc.zip') as myzip:
      with myzip.open(isolate + '_L001_' + read + '_001_corrected_fastqc/summary.txt') as myfile:
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

  print(isolate, SPECIES, PCT, READS, ST_SCHEME, ST, CONTIGS, N50, LARGEST, SIZE, COVERAGE, WARN_PRE_FINAL, FAIL_PRE_FINAL, WARN_POST_FINAL, FAIL_POST_FINAL, sep = ',')
