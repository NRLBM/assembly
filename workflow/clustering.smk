import sys
import os

configfile: "workflow/config/config.yaml"

def get_samples_from_qc_report(qc_report_file, timestamp):
  tmp_list = []
  full_path_qc_file = ''.join(['backup/', timestamp, '/', qc_report_file])
  with open(full_path_qc_file) as qc_file:
    lines = qc_file.readlines()[1:]
  for line in lines:
    tmp_list.append(line.split('\t')[0])
  return tmp_list

ECOLI_SAMPLES = get_samples_from_qc_report('qc_report_ecoli.tsv', config["timestamp"])
NMEN_SAMPLES = get_samples_from_qc_report('qc_report_nmen.tsv', config["timestamp"])
SPYO_SAMPLES = get_samples_from_qc_report('qc_report_spyo.tsv', config["timestamp"])

CLUSTERS_NMEN, = glob_wildcards('/'.join(["tmp_data", config["timestamp"], "cluster_lists/{cluster_ref}.txt"]))

#print(CLUSTERS_NMEN)

all_output = []

if len(CLUSTERS_NMEN) > 0:
#  all_output.append(expand("tmp_data/{timestamp}/microreact/{cluster_ref}/{cluster_ref}.nwk", cluster_ref=CLUSTERS_NMEN, timestamp=config["timestamp"]))
#  all_output.append(expand("tmp_data/{timestamp}/microreact/{cluster_ref}/{cluster_ref}.dot", cluster_ref=CLUSTERS_NMEN, timestamp=config["timestamp"]))
#  all_output.append(expand("tmp_data/{timestamp}/microreact/{cluster_ref}/{cluster_ref}.csv", cluster_ref=CLUSTERS_NMEN, timestamp=config["timestamp"]))
  all_output.append(expand("backup/{timestamp}/microreact/{cluster_ref}", cluster_ref=CLUSTERS_NMEN, timestamp=config["timestamp"])),
  all_output.append(expand("workflow/db/mash_databases/{timestamp}.msh", timestamp=config["timestamp"]))

if len(all_output) == 0:
  print("No cluster lists could be found so no clustering analysis is performed.")
  sys.exit(0)

rule all:
  input:
    all_output

rule ska_fasta_nmen:
  input:
    "tmp_data/{timestamp}/cluster_lists/{cluster_ref}.txt"
  output:
    directory("tmp_data/{timestamp}/ska_fasta/{cluster_ref}")
  conda:
    "envs/ska.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/ska_fasta_Nmen/{cluster_ref}.log"
  threads: 8
  shell:
    """
    bash workflow/scripts/ska_wrapper.sh {input} {output} 2>&1>{log}
    """

rule ska_align_nmen:
  input:
    "tmp_data/{timestamp}/ska_fasta/{cluster_ref}"
  output:
    "tmp_data/{timestamp}/ska_align/{cluster_ref}_variants.aln"
  conda:
    "envs/ska.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/ska_align_Nmen/{cluster_ref}.log"
  params:
    proportion = config["ska"]["proportion"],
    outfile = "tmp_data/{timestamp}/ska_align/{cluster_ref}"
  threads: 16
  shell:
    """
    ska align -v -o {params.outfile} -p {params.proportion} {input}/* 2>&1>{log}
    """

rule fasttree_nmen:
  input:
    "tmp_data/{timestamp}/ska_align/{cluster_ref}_variants.aln"
  output:
    "tmp_data/{timestamp}/microreact/{cluster_ref}/{cluster_ref}.nwk"
  conda:
    "envs/fasttree.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/fasttree_Nmen/{cluster_ref}.log"
  threads: 16
  shell:
    """
    fasttree -nt {input} > {output}
    """

rule snp_dists_nmen:
  input:
    "tmp_data/{timestamp}/ska_align/{cluster_ref}_variants.aln"
  output:
    "tmp_data/{timestamp}/snp_dists_out/{cluster_ref}.tsv"
  conda:
    "envs/snpdists.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/snp_dists_Nmen/{cluster_ref}.log"
  threads: 16
  shell:
    """
    snp-dists -m {input} > {output}
    """

rule convert_snps_to_network_nmen:
  input:
    "tmp_data/{timestamp}/snp_dists_out/{cluster_ref}.tsv"
  output:
    "tmp_data/{timestamp}/microreact/{cluster_ref}/{cluster_ref}.dot"
  conda:
    "envs/convert_snps_to_network.yaml"
  params:
    threshold = config["convert_snps_to_network"]["threshold"]
  log:
    "slurm/snakemake_logs/{timestamp}/convert_snps_to_network_Nmen/{cluster_ref}.log"
  shell:
    """
    python workflow/scripts/convert_SNPs_to_network.py -t {params.threshold} -i {input} -o {output}
    """

rule mark_new_samples_nmen:
  input:
    "tmp_data/{timestamp}/snp_dists_out/{cluster_ref}.tsv"
  output:
    "tmp_data/{timestamp}/microreact/{cluster_ref}/{cluster_ref}.csv"
  conda:
    "envs/convert_snps_to_network.yaml"
  params:
    samples = NMEN_SAMPLES
  log:
    "slurm/snakemake_logs/{timestamp}/mark_new_samples_nmen/{cluster_ref}.log"
  shell:
    """
    python workflow/scripts/mark_new_samples.py -i {input} -o {output} {params.samples}
    """

rule backup_microreact_data_nmen:
  input:
    csv = "tmp_data/{timestamp}/microreact/{cluster_ref}/{cluster_ref}.csv",
    dot = "tmp_data/{timestamp}/microreact/{cluster_ref}/{cluster_ref}.dot",
    nwk = "tmp_data/{timestamp}/microreact/{cluster_ref}/{cluster_ref}.nwk"
  output:
    directory("backup/{timestamp}/microreact/{cluster_ref}")
  shell:
    """
    mkdir -p {output}
    cp {input} {output}
    """

rule add_mash_db_nmen:
  input:
    genomes = expand("output/{timestamp}/genomes/{sample}.fasta", timestamp=config["timestamp"], sample=NMEN_SAMPLES),
    dummy = expand("backup/{timestamp}/microreact/{cluster_ref}", timestamp=config["timestamp"], cluster_ref=CLUSTERS_NMEN),
  output:
    expand("workflow/db/mash_databases/{timestamp}.msh", timestamp=config["timestamp"])
  conda:
    "envs/mash.yaml"
  log:
    expand("slurm/snakemake_logs/{timestamp}/add_mash_db_nmen.log", timestamp=config["timestamp"])
  params:
    prefix = "workflow/db/mash_databases",
    timestamp = config["timestamp"]
  threads: 16
  shell:
    """
    mash sketch -o {params.prefix}/{params.timestamp} {input.genomes}
    """

