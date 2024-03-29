#!/usr/bin/env python3

__author__ = "Boas van der Putten"
__version__ = "0.0.1"

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Summarise typed isolates')

parser.add_argument('samples', nargs='+', help="List of sample names to summarise", type=str)
parser.add_argument('--species', dest='species', help="Species to summarise", choices=['Escherichia coli', 'Neisseria meningitidis', 'Ecoli', 'Nmen', 'Streptococcus pyogenes', 'Spyo'], type=str, required=True)
parser.add_argument('--qc', dest='qc', help="Path to QC report", type=str, required=True)
parser.add_argument('--mlst', dest="mlst", help="mlst output directory", type=str)
parser.add_argument('--gmats', dest="gmats", help="gMATS output directory", type=str, default="gMATS_Nmen")
parser.add_argument('--meningotype', dest="meningotype", help="meningotype output directory", type=str, default="meningotype_Nmen")
parser.add_argument('--amrfinder', dest="amrfinder", help="AMRfinder output directory", type=str)
parser.add_argument('--ectyper', dest="ectyper", help="ECtyper output directory", type=str, default="ectyper_Ecoli")
parser.add_argument('--fimtyper', dest="fimtyper", help="ABRicate Fimtyper output directory", type=str, default="ABRicate_fimH_Ecoli")
parser.add_argument('--vfdb', dest="vfdb", help="ABRicate VFDB output directory", type=str)
parser.add_argument('--emmtyper', dest="emmtyper", help="emmtyper output directory", type=str, default="emmtyper_Spyo")
parser.add_argument('--M1UK', dest="M1UK", help="M1UK mapping output directory", type=str, default="M1UK_Spyo")
parser.add_argument('--timestamp', dest="timestamp", help="Timestamp of analysis", type=str, required=True)
parser.add_argument('--output', dest='output', help="Output file", type=str, required=True)
parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))

args = parser.parse_args()

if (args.species == "Escherichia coli") or (args.species == "Ecoli"):
  args.species = "Ecoli"
  assert args.mlst is not None
  assert args.amrfinder is not None
  assert args.ectyper is not None
  assert args.fimtyper is not None
  assert args.vfdb is not None
elif (args.species == "Neisseria meningitidis") or (args.species == "Nmen"):
  args.species = "Nmen"
  assert args.mlst is not None
  assert args.amrfinder is not None
  assert args.vfdb is not None
  assert args.gmats is not None
  assert args.meningotype is not None
elif (args.species == "Streptococcus pyogenes") or (args.species == "Spyo"):
  args.species = "Spyo"
  assert args.mlst is not None
  assert args.emmtyper is not None
  assert args.M1UK is not None

# Store samples as isolate_list
isolate_list = list(args.samples)

# Make output prefix based on timestamp
output_prefix = 'tmp_data/' + str(args.timestamp) + '/'

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

def parse_amrfinder(gene, lines, WT):
  tmp_list = []
  for line in lines:
    if gene in line.split('\t')[5]:
      tmp_list.append(line.split('\t')[5])
  if len(tmp_list) == 0:
    tmp_list.append(WT)
  output_string = '|'.join(tmp_list)
  return output_string

def parse_vfdb(gene, lines):
  tmp_list = []
  for line in lines:
    if gene in line.split('\t')[5]:
      tmp_list.append(line.split('\t')[5])
  if len(tmp_list) == 0:
    output_string = 'None'
  else:
    tmp_list.sort()
    output_string = gene
  for i in tmp_list:
    output_string = output_string + str(i[-1:])
  return output_string

def parse_M1UK(row):
  if row['emm_type'] != 'EMM1.0':
    return 'non-M1'
  elif row['exact_match'] == 27:
    return 'M1UK'
  elif row['exact_match'] > 14:
    return 'check_manually'
  else:
    return 'M1global'

def add_qc_results(full_df, qc_df_path):
  qc_df = pd.read_csv(qc_df_path, sep='\t')
  full_df = qc_df.merge(full_df, on='Isolate', how='right', validate='one_to_one')
  full_df = full_df.fillna('-')
  return full_df

def process_ecoli(args):
  full_df = pd.DataFrame()

  tuple_resistance_genes = ('blaCTX-M', 'blaOXA', 'blaCMY', 'gyrA_S83', 'gyrA_D87', 'parC_S80')
  tuple_VFDB_genes = ('asl', 'cnf', 'fde', 'ibe', 'kps', 'neu', 'sfa', 'hly', 'iro', 'iut', 'iuc', 'pap', 'ybt')

  for isolate in isolate_list:
    # Get ST
    filepath = output_prefix + args.mlst + '/' + isolate + ".tsv"
    isolate_df = pd.read_csv(filepath, sep='\t')

    # Get AMR genes and put lines in a list to reuse
    lines = customreadlines(output_prefix + args.amrfinder + '/' + isolate + ".tsv")

    for resistance_gene in tuple_resistance_genes:
      isolate_df[resistance_gene] = parse_amrfinder(resistance_gene, lines, '-')

    ############## ECTYPER ################
    filepath = output_prefix + args.ectyper + '/' + isolate + "/output.tsv"
    ectyper_df = pd.read_csv(filepath, sep='\t')
    ectyper_df = ectyper_df.drop(['Name', 'Species', 'QC', 'AlleleKeys', 'GeneContigNames', 'GeneRanges', 'GeneLengths'], axis=1)

    isolate_df = pd.concat([isolate_df, ectyper_df], axis=1)

    # ABRicate fimH typing
    lines = customreadlines(output_prefix + args.fimtyper + '/' + isolate + '.tsv')
    isolate_df['fimH'] = parse_amrfinder('fimH', lines, 'None')

    # ABRicate VFDB typing
    lines = customreadlines(output_prefix + args.vfdb + "/" + isolate + ".tsv")

    for VFDB_gene in tuple_VFDB_genes:
      isolate_df[VFDB_gene] = parse_vfdb(VFDB_gene, lines)

    full_df = full_df.append(isolate_df)

  # Read QC table and merge
  full_df = add_qc_results(full_df, args.qc)

  return full_df

def process_nmen(args):
  full_df = pd.DataFrame()

  tuple_resistance_genes = ('penA_A311', 'penA_I312', 'penA_V316', 'penA_D346D', 'penA_T483', 'penA_A501', 'penA_F504', 'penA_A510', 'penA_N512', 'penA_A516', 'penA_G542', 'penA_G545', 'penA_P551', 'folP_R228', 'mtrR_A39', 'mtrR_G45', 'gyrA_S91', 'gyrA_D95', 'gyrB_D429', 'gyrB_K450', 'parC_D86', 'parC_S87', 'parC_S88', 'parC_E91', 'parE_G410', 'ponA_L421', 'porB1b_G120', 'porB1b_A121', 'rplD_G70', 'rplV_K83KKGPSL', 'rplV_K90KARA', 'rpoB_G158', 'rpoB_P157', 'rpoB_R201', 'rpoD_DDDA92de', 'rpoD_E98', 'rpsE_K28', 'rpsE_T24', 'rpsE_V27de', 'rpsJ_V57', '16S_C1191', '23S_A2045', '23S_A2057', '23S_A2145', '23S_C2597', 'macAB_G-48', 'mtrC_C-120', 'mtrR_A-53', 'mtrR_A-56', 'mtrR_G-131', 'norM_A-7', 'norM_C-104')
  tuple_VFDB_genes = ('ctr', 'far', 'fbp', 'fHbp', 'hmb', 'hpu', 'iga', 'kat', 'lbp', 'lip', 'mnt', 'mtr', 'opc', 'pil', 'por', 'sia')

  for isolate in isolate_list:
    # Get ST
    filepath = output_prefix + args.mlst + '/' + isolate + ".tsv"
    isolate_df = pd.read_csv(filepath, sep='\t')

    ### GMATS
    filepath = output_prefix + args.gmats + '/' + isolate + ".tsv"
    gmats_df = pd.read_csv(filepath, sep='\t')
    gmats_df = gmats_df.drop(['Isolate'], axis=1)

    isolate_df = pd.concat([isolate_df, gmats_df], axis=1)

    ### meningotype
    filepath = output_prefix + args.meningotype + '/' + isolate + ".tsv"
    meningotype_df = pd.read_csv(filepath, sep='\t')
    isolate_df['Serogroup'] = meningotype_df.at[0,'SEROGROUP']
    isolate_df['porA'] = meningotype_df.at[0,'PorA']
    isolate_df['fetA'] = meningotype_df.at[0,'FetA']

    # Get AMR genes and put lines in a list to reuse
    lines = customreadlines(output_prefix + args.amrfinder + '/' + isolate + ".tsv")

    for resistance_gene in tuple_resistance_genes:
      isolate_df[resistance_gene] = parse_amrfinder(resistance_gene, lines, '-')

    # ABRicate VFDB typing
    lines = customreadlines(output_prefix + args.vfdb + "/" + isolate + ".tsv")

    for VFDB_gene in tuple_VFDB_genes:
      isolate_df[VFDB_gene] = parse_vfdb(VFDB_gene, lines)

    full_df = full_df.append(isolate_df)

  full_df = add_qc_results(full_df, args.qc)

  return full_df

def process_spyo(args):
  full_df = pd.DataFrame()

  for isolate in isolate_list:
    # Get ST
    filepath = output_prefix + args.mlst + '/' + isolate + ".tsv"
    isolate_df = pd.read_csv(filepath, sep='\t')

    # Get emmtyper results
    filepath_emmtyper = output_prefix + args.emmtyper + '/' + isolate + '.tsv'
    emmtyper_df = pd.read_csv(filepath_emmtyper, sep='\t', names=['Isolate_raw', 'nr_blast_hits', 'nr_clusters', 'emm_type', 'emm_gene_position', 'possible_emm_types', 'possible_emm_gene_locations', 'emm_cluster'])
    emmtyper_df['Isolate'] = os.path.splitext(os.path.basename(emmtyper_df['Isolate_raw'][0]))[0]
    emmtyper_df = emmtyper_df.drop(['Isolate_raw', 'nr_blast_hits', 'nr_clusters', 'possible_emm_types', 'possible_emm_gene_locations'], axis=1)

    isolate_df = isolate_df.merge(emmtyper_df, on='Isolate', how='left')

    ###### M1UK
    filepath_M1UK = output_prefix + args.M1UK + '/' + isolate + '.tsv'
    M1UK_df = pd.read_csv(filepath_M1UK, sep='\t')
    M1UK_df = M1UK_df.drop(['115646_C', '116162_A', '116163_C', '250832_T', '513254_G', '528360_A', '563631_G', '613633_T', '626494_G', '661707_G', '730823_C', '784467_T', '819098_G', '923079_G', '942633_G', '983438_G', '1082253_C', '1238124_G', '1238673_G', '1251193_G', '1373176_C', '1407497_C', '1446116_C', '1535209_A', '1702540_C', '1734749_G', '1828734_G'], axis=1)

    ### Based on emmtyper results and M1UK SNPs, decide non-M1/M1global/M1UK
    isolate_df = isolate_df.merge(M1UK_df, on='Isolate', how='left')

    isolate_df['M1_type'] = isolate_df.apply(parse_M1UK, axis=1)

    full_df = full_df.append(isolate_df)

  full_df = add_qc_results(full_df, args.qc)

  return full_df


if __name__ == "__main__":
  if args.species == "Ecoli":
    full_df = process_ecoli(args)
  elif args.species == "Nmen":
    full_df = process_nmen(args)
  elif args.species == "Spyo":
    full_df = process_spyo(args)

  full_df.to_csv(args.output, sep = '\t', index=False)

