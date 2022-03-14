#!/usr/bin/env python3

def get_sample_name(lines):
  for line in lines:
    if line.split('\t')[0] == '#CHROM':
      return line.split('\t')[9].rstrip('\n')

def clean_lines(lines):
  out = list()
  for line in lines:
    if line[0] == '#':
      continue
    line_cols = line.rstrip('\n').split('\t')
    out.append(line_cols)
  return out

def get_ref_dict(lines):
  ref_dict = {}
  for line in lines:
    ref_dict[line[1]] = (line[3], line[4])
  return ref_dict

def compare_sample_ref(lines, ref_dict):
  exact_sample_count_match = 0
  variant_sample_count_match = 0
  sample_count_nonmatch = 0
  SNP_identities = []
  for position, identities in ref_dict.items():
    non_match_increment = 1
    SNP_identity = identities[0]
    for line in lines:
      if line[1] == position:
        non_match_increment = 0
        SNP_ref = identities[0]
        assert SNP_ref == line[3], f"SNP REF of position {position} should be {SNP_ref} but is {line[3]}. Is NC_007297.2 the reference?"
        SNP_identity = line[4]
        if SNP_identity == identities[1]:
          exact_sample_count_match += 1
        else:
          variant_sample_count_match += 1
    sample_count_nonmatch += non_match_increment
    SNP_identities.append(SNP_identity)
  SNP_identities = '\t'.join(SNP_identities)
  return exact_sample_count_match, variant_sample_count_match, sample_count_nonmatch, SNP_identities

def get_SNP_ref_identies(dict):
  list_colnames = []
  for position, identities in dict.items():
    colname = '_'.join([position, identities[0]])
    list_colnames.append(colname)
  return list_colnames

def print_output(args, sample_name, exact_match, variant_match, no_match, SNP_identities, ref_dict):
  if args.extended:
    args.header = True
    outstring =  '\t'.join([str(x) for x in [sample_name, exact_match, variant_match, no_match, SNP_identities]])
  else:
    outstring =  '\t'.join([str(x) for x in [sample_name, exact_match, variant_match, no_match]])
  
  if args.header:
    if args.extended:
      SNP_ref_identities = get_SNP_ref_identies(ref_dict)
      headerstring = '\t'.join(['Isolate', 'exact_match', 'variant_match', 'non_match'] + SNP_ref_identities)
    else:
      headerstring = '\t'.join(['Isolate', 'exact_match', 'variant_match', 'non_match'])
    print(headerstring)
  
  print(outstring)

def main(args):
  if args.input == sys.stdin:
    sample_lines = sys.stdin.readlines()
  else:
    with open(args.input) as input_file:
      sample_lines = input_file.readlines()
  
  with open(args.ref) as ref_file:
    ref_lines = ref_file.readlines()
  
  sample_name = get_sample_name(sample_lines)
  ref_lines_clean = clean_lines(ref_lines)
  sample_lines_clean = clean_lines(sample_lines)
  ref_dict = get_ref_dict(ref_lines_clean)
  exact_match, variant_match, no_match, SNP_identities = compare_sample_ref(sample_lines_clean, ref_dict)
  print_output(args, sample_name, exact_match, variant_match, no_match, SNP_identities, ref_dict)

if __name__ == '__main__':
  import argparse
  import sys

  # Parse argument
  parser = argparse.ArgumentParser(description='Count lineage-specific M1UK SNPs.')
  parser.add_argument('-i', '--input', help="Input VCF file, filtered for M1UK-specific SNPs and decomposed using vt decompose_blocksub", dest='input', default=sys.stdin)
  parser.add_argument('-r', '--ref', help="Reference VCF file containing lineage-specific SNPs", dest='ref', required=True, type=str)
  parser.add_argument('--extended', help="Print SNP identities", dest='extended', action='store_true')
  parser.add_argument('--header', help="Print header", dest='header', action='store_true')
  parser.add_argument('-o', '--output', help="Output file", dest='output', default=sys.stderr)
  args = parser.parse_args()

  main(args)
