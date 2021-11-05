import sys

configfile: "workflow/config/config.yaml"

def get_samples_from_qc_report(qc_report_file, timestamp):
  tmp_list = []
  full_path_qc_file = ''.join(['backup/', timestamp, '/', qc_report_file])
  with open(full_path_qc_file) as qc_file:
    lines = qc_file.readlines()
  for line in lines:
    if line.split('\t')[-1].rstrip('\n') == 'PASS':
      tmp_list.append(line.split('\t')[0])
  return tmp_list

ECOLI_SAMPLES = get_samples_from_qc_report('qc_report_ecoli.tsv', config["timestamp"])
NMEN_SAMPLES = get_samples_from_qc_report('qc_report_nmen.tsv', config["timestamp"])

all_output = []

if len(ECOLI_SAMPLES) > 0:
#  all_output.append(expand("tmp_data/{timestamp}/ectyper_Ecoli/{sample}", sample=ECOLI_SAMPLES, timestamp=config["timestamp"]))
#  all_output.append(expand("tmp_data/{timestamp}/MLST_Ecoli/{sample}.tsv", sample=ECOLI_SAMPLES, timestamp=config["timestamp"]))
#  all_output.append(expand("tmp_data/{timestamp}/AMRfinder_Ecoli/{sample}.tsv", sample=ECOLI_SAMPLES, timestamp=config["timestamp"]))
#  all_output.append(expand("tmp_data/{timestamp}/ABRicate_fimH_Ecoli/{sample}.tsv", sample=ECOLI_SAMPLES, timestamp=config["timestamp"]))
#  all_output.append(expand("tmp_data/{timestamp}/ABRicate_VFDB_Ecoli/{sample}.tsv", sample=ECOLI_SAMPLES, timestamp=config["timestamp"]))
  all_output.append(expand("output/{timestamp}/Ecoli_typing_summary.xlsx", timestamp=config["timestamp"]))
  all_output.append(expand("backup/{timestamp}/Ecoli_typing_summary.tsv", timestamp=config["timestamp"]))
  all_output.append(expand("backup/{timestamp}/Ecoli_typing/{sample}_typing.tar.gz", sample=ECOLI_SAMPLES, timestamp=config["timestamp"]))

if len(NMEN_SAMPLES) > 0:
#  all_output.append(expand("tmp_data/{timestamp}/MLST_Nmen/{sample}.tsv", sample=NMEN_SAMPLES, timestamp=config["timestamp"]))
#  all_output.append(expand("tmp_data/{timestamp}/gMATS_Nmen/{sample}.json", sample=NMEN_SAMPLES, timestamp=config["timestamp"]))
#  all_output.append(expand("tmp_data/{timestamp}/gMATS_Nmen/{sample}.tsv", sample=NMEN_SAMPLES, timestamp=config["timestamp"]))
#  all_output.append(expand("tmp_data/{timestamp}/AMRfinder_Nmen/{sample}.tsv", sample=NMEN_SAMPLES, timestamp=config["timestamp"]))
#  all_output.append(expand("tmp_data/{timestamp}/ABRicate_VFDB_Nmen/{sample}.tsv", sample=NMEN_SAMPLES, timestamp=config["timestamp"]))
  all_output.append(expand("output/{timestamp}/Nmen_typing_summary.xlsx", timestamp=config["timestamp"]))
  all_output.append(expand("backup/{timestamp}/Nmen_typing_summary.tsv", timestamp=config["timestamp"]))
  all_output.append(expand("backup/{timestamp}/Nmen_typing/{sample}_typing.tar.gz", sample=NMEN_SAMPLES, timestamp=config["timestamp"]))

if len(all_output) == 0:
  print("No samples to process", file=sys.stderr)

rule all:
  input:
    all_output

# Escherichia coli typing
rule MLST_ecoli:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta"
  output:
    json = temp("tmp_data/{timestamp}/MLST_Ecoli/{sample}.json")
  threads: 8
  params:
    scheme = config['pubmlst']['ecoli_mlst_scheme']
  conda:
    "envs/python.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/MLST_Ecoli_PubMLST_API/{sample}.log"
  shell:
    """
    python workflow/scripts/type_pubmlst_api.py --genome {input.genome} --api-url {params.scheme} 1> {output.json} 2>{log}
    """

rule convert_MLST_ecoli:
  input:
    "tmp_data/{timestamp}/MLST_Ecoli/{sample}.json"
  output:
    temp("tmp_data/{timestamp}/MLST_Ecoli/{sample}.tsv")
  threads: 1
  conda:
    "envs/python.yaml"
  params:
    loci = config['pubmlst']['ecoli_mlst_loci']
  log:
    "slurm/snakemake_logs/{timestamp}/MLST_Ecoli_convert/{sample}.log"
  shell:
    """
    python workflow/scripts/convert_MLST.py --input {input} --output {output} --loci {params.loci}
    """

rule AMRFinder_ecoli:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta"
  output: 
    temp("tmp_data/{timestamp}/AMRfinder_Ecoli/{sample}.tsv")
  conda:
    "envs/amrfinder.yaml"
  params: 
    organism = "Escherichia"
  threads: 16
  log:
    "slurm/snakemake_logs/{timestamp}/AMRfinder_Ecoli/{sample}.log"
  shell:
    """
    amrfinder -u
    amrfinder --threads {threads} --nucleotide {input} --organism {params.organism} --output {output} 2>&1>{log}    
    """

rule ECtyper_ecoli:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta"
  output:
    directory("tmp_data/{timestamp}/ectyper_Ecoli/{sample}")
  conda:
    "envs/ectyper.yaml"
  params:
    general = config["ectyper"]["general"]
  threads: 8
  log:
    "slurm/snakemake_logs/{timestamp}/ECtyper/{sample}.log"
  shell:
    """
    ectyper -i {input.genome} -c {threads} -o {output} 2>&1>{log}
    """

rule ABRicate_fimH_Ecoli:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta"
  output:
    tsv = "tmp_data/{timestamp}/ABRicate_fimH_Ecoli/{sample}.tsv"
  conda:
    "envs/abricate.yaml"
  params:
    datadir = config["abricate"]["fimH_datadir"],
    db = config["abricate"]["fimH_db"],
    minid = config["abricate"]["minid"],
    mincov = config["abricate"]["mincov"]
  threads: 1
  log:
    "slurm/snakemake_logs/{timestamp}/ABRicate_fimH_Ecoli/{sample}.log"
  shell:
    """
    abricate --db {params.db} --datadir {params.datadir} --threads {threads} --nopath --minid {params.minid} --mincov {params.mincov} {input.genome} 1>{output.tsv} 2>{log}
    """

rule ABRicate_VFDB_Ecoli:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta"
  output:
    tsv = "tmp_data/{timestamp}/ABRicate_VFDB_Ecoli/{sample}.tsv"
  conda:
    "envs/abricate.yaml"
  params:
    db = config["abricate"]["VFDB_db"],
    minid = config["abricate"]["minid"],
    mincov = config["abricate"]["mincov"]
  threads: 1
  log:
    "slurm/snakemake_logs/{timestamp}/ABRicate_VFDB_Ecoli/{sample}.log"
  shell:
    """
    abricate --db {params.db} --threads {threads} --nopath --minid {params.minid} --mincov {params.mincov} {input.genome} 1>{output.tsv} 2>{log}
    """

# Find all files to backup and create a tarball with a timestamp
rule backup_data_ecoli:
  input:
    ectyper = "tmp_data/{timestamp}/ectyper_Ecoli/{sample}",
    MLST_tsv = "tmp_data/{timestamp}/MLST_Ecoli/{sample}.tsv",
    MLST_json = "tmp_data/{timestamp}/MLST_Ecoli/{sample}.json",
    AMRfinder = "tmp_data/{timestamp}/AMRfinder_Ecoli/{sample}.tsv",
    ABRicate_fimH = "tmp_data/{timestamp}/ABRicate_fimH_Ecoli/{sample}.tsv",
    ABRicate_VFDB = "tmp_data/{timestamp}/ABRicate_VFDB_Ecoli/{sample}.tsv",
  output:
    "backup/{timestamp}/Ecoli_typing/{sample}_typing.tar.gz"
  params:
    sample = "{sample}",
    timestamp = config["timestamp"]
  threads: 16
  shell:
    """
    mkdir -p tmp_data/.backup/{params.sample}/MLST_PubMLST
    mkdir -p tmp_data/.backup/{params.sample}/AMRfinder
    mkdir -p tmp_data/.backup/{params.sample}/ABRicate_fimH
    mkdir -p tmp_data/.backup/{params.sample}/ABRicate_VFDB
    cp {input.MLST_tsv} tmp_data/.backup/{params.sample}/MLST_PubMLST
    cp {input.MLST_json} tmp_data/.backup/{params.sample}/MLST_PubMLST
    cp {input.AMRfinder} tmp_data/.backup/{params.sample}/AMRfinder
    cp {input.ABRicate_fimH} tmp_data/.backup/{params.sample}/ABRicate_fimH
    cp {input.ABRicate_VFDB} tmp_data/.backup/{params.sample}/ABRicate_VFDB
    cp -r {input.ectyper} tmp_data/.backup/{params.sample}/ectyper
    tar zcvf {output} -C tmp_data/.backup {params.sample}
    rm -rf tmp_data/.backup
    """

rule typing_summary_csv_Ecoli:
  input:
    expand("tmp_data/{timestamp}/ectyper_Ecoli/{sample}", sample=ECOLI_SAMPLES, timestamp=config["timestamp"]),
    expand("tmp_data/{timestamp}/MLST_Ecoli/{sample}.tsv", sample=ECOLI_SAMPLES, timestamp=config["timestamp"]),
    expand("tmp_data/{timestamp}/AMRfinder_Ecoli/{sample}.tsv", sample=ECOLI_SAMPLES, timestamp=config["timestamp"]),
    expand("tmp_data/{timestamp}/ABRicate_fimH_Ecoli/{sample}.tsv", sample=ECOLI_SAMPLES, timestamp=config["timestamp"]),
    expand("tmp_data/{timestamp}/ABRicate_VFDB_Ecoli/{sample}.tsv", sample=ECOLI_SAMPLES, timestamp=config["timestamp"])
  output:
    "backup/{timestamp}/Ecoli_typing_summary.tsv"
  conda:
    "envs/python.yaml"
  params:
    timestamp = config["timestamp"],
    samples = ECOLI_SAMPLES
  threads: 16
  log:
    "slurm/snakemake_logs/{timestamp}/Ecoli_typing_summary.log"
  shell:
    """
    python workflow/scripts/typing_summary.py --species Ecoli --timestamp {params.timestamp} --output {output} --mlst MLST_Ecoli --amrfinder AMRfinder_Ecoli --vfdb ABRicate_VFDB_Ecoli {params.samples} 2>&1>{log}
    """

# Convert the summary file (tsv format) to an Excel file and save this version in the output folder
rule typing_summary_to_xlsx_Ecoli:
  input:
    "backup/{timestamp}/Ecoli_typing_summary.tsv"
  output:
    "output/{timestamp}/Ecoli_typing_summary.xlsx"
  conda:
    "envs/python.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/Ecoli_typing_summary_to_xlsx.log"
  threads: 16
  shell:
    """
    python workflow/scripts/summary_to_xlsx.py -d '\t' {input} {output} 2>&1>{log}
    """

# Neisseria meningitidis typing
rule MLST_nmen:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta"
  output:
    json = temp("tmp_data/{timestamp}/MLST_Nmen/{sample}.json")
  threads: 8
  params:
    scheme = config['pubmlst']['nmen_mlst_scheme']
  conda:
    "envs/python.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/MLST_Nmen_PubMLST_API/{sample}.log"
  shell:
    """
    python workflow/scripts/type_pubmlst_api.py --genome {input.genome} --api-url {params.scheme} 1> {output.json} 2>{log}
    """

rule convert_MLST_nmen:
  input:
    "tmp_data/{timestamp}/MLST_Nmen/{sample}.json"
  output:
    temp("tmp_data/{timestamp}/MLST_Nmen/{sample}.tsv")
  threads: 1
  conda:
    "envs/python.yaml"
  params:
    loci = config['pubmlst']['nmen_mlst_loci']
  log:
    "slurm/snakemake_logs/{timestamp}/MLST_Nmen_convert/{sample}.log"
  shell:
    """
    python workflow/scripts/convert_MLST.py --input {input} --output {output} --loci {params.loci}
    """

rule AMRFinder_nmen:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta"
  output:
    temp("tmp_data/{timestamp}/AMRfinder_Nmen/{sample}.tsv")
  conda:
    "envs/amrfinder.yaml"
  params:
    organism = "Neisseria"
  threads: 16
  log:
    "slurm/snakemake_logs/{timestamp}/AMRfinder_Nmen/{sample}.log"
  shell:
    """
    amrfinder -u
    amrfinder --threads {threads} --nucleotide {input} --organism {params.organism} --output {output} 2>&1>{log}
    """

rule gMATS_nmen:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta"
  output:
    json = temp("tmp_data/{timestamp}/gMATS_Nmen/{sample}.json")
  threads: 8
  params:
    scheme = config['pubmlst']['nmen_vaccine_antigen_scheme']
  conda:
    "envs/python.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/gMATS_Nmen_PubMLST_API/{sample}.log"
  shell:
    """
    python workflow/scripts/type_pubmlst_api.py --genome {input.genome} --api-url {params.scheme} 1> {output.json} 2>{log}
    """

rule gMATS_coverage_Nmen:
  input:
    json = "tmp_data/{timestamp}/gMATS_Nmen/{sample}.json"
  output:
    tsv = "tmp_data/{timestamp}/gMATS_Nmen/{sample}.tsv"
  conda:
    "envs/python.yaml"
  threads: 1
  log:
    "slurm/snakemake_logs/{timestamp}/gMATS_coverage_Nmen/{sample}.log"
  shell:
    """
    python workflow/scripts/gMATS_coverage.py --input {input.json} --output {output.tsv} 2>&1>{log}
    """

rule ABRicate_VFDB_Nmen:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta"
  output:
    tsv = "tmp_data/{timestamp}/ABRicate_VFDB_Nmen/{sample}.tsv"
  conda:
    "envs/abricate.yaml"
  params:
    db = config["abricate"]["VFDB_db"],
    minid = config["abricate"]["minid"],
    mincov = config["abricate"]["mincov"]
  threads: 1
  log:
    "slurm/snakemake_logs/{timestamp}/ABRicate_VFDB_Ecoli/{sample}.log"
  shell:
    """
    abricate --db {params.db} --threads {threads} --nopath --minid {params.minid} --mincov {params.mincov} {input.genome} 1>{output.tsv} 2>{log}
    """

# Find all files to backup and create a tarball with a timestamp
rule backup_data_nmen:
  input:
    gMATS_tsv = "tmp_data/{timestamp}/gMATS_Nmen/{sample}.tsv",
    gMATS_json = "tmp_data/{timestamp}/gMATS_Nmen/{sample}.json",
    MLST_tsv = "tmp_data/{timestamp}/MLST_Nmen/{sample}.tsv",
    MLST_json = "tmp_data/{timestamp}/MLST_Nmen/{sample}.json",
    AMRfinder = "tmp_data/{timestamp}/AMRfinder_Nmen/{sample}.tsv",
    ABRicate_VFDB = "tmp_data/{timestamp}/ABRicate_VFDB_Nmen/{sample}.tsv",
  output:
    "backup/{timestamp}/Nmen_typing/{sample}_typing.tar.gz"
  params:
    sample = "{sample}",
    timestamp = config["timestamp"]
  threads: 16
  shell:
    """
    mkdir -p tmp_data/.backup/{params.sample}/MLST_PubMLST
    mkdir -p tmp_data/.backup/{params.sample}/AMRfinder
    mkdir -p tmp_data/.backup/{params.sample}/gMATS
    mkdir -p tmp_data/.backup/{params.sample}/ABRicate_VFDB
    cp {input.gMATS_tsv} tmp_data/.backup/{params.sample}/gMATS
    cp {input.gMATS_json} tmp_data/.backup/{params.sample}/gMATS
    cp {input.MLST_tsv} tmp_data/.backup/{params.sample}/MLST_PubMLST
    cp {input.MLST_json} tmp_data/.backup/{params.sample}/MLST_PubMLST
    cp {input.AMRfinder} tmp_data/.backup/{params.sample}/AMRfinder
    cp {input.ABRicate_VFDB} tmp_data/.backup/{params.sample}/ABRicate_VFDB
    tar zcvf {output} -C tmp_data/.backup {params.sample}
    rm -rf tmp_data/.backup
    """

rule typing_summary_csv_Nmen:
  input:
    expand("tmp_data/{timestamp}/gMATS_Nmen/{sample}.tsv", sample=NMEN_SAMPLES, timestamp=config["timestamp"]),
    expand("tmp_data/{timestamp}/MLST_Nmen/{sample}.tsv", sample=NMEN_SAMPLES, timestamp=config["timestamp"]),
    expand("tmp_data/{timestamp}/AMRfinder_Nmen/{sample}.tsv", sample=NMEN_SAMPLES, timestamp=config["timestamp"]),
    expand("tmp_data/{timestamp}/ABRicate_VFDB_Nmen/{sample}.tsv", sample=NMEN_SAMPLES, timestamp=config["timestamp"])
  output:
    "backup/{timestamp}/Nmen_typing_summary.tsv"
  conda:
    "envs/python.yaml"
  params:
    timestamp = config["timestamp"],
    samples = NMEN_SAMPLES
  threads: 16
  log:
    "slurm/snakemake_logs/{timestamp}/Nmen_typing_summary.log"
  shell:
    """
    python workflow/scripts/typing_summary.py --species Nmen --timestamp {params.timestamp} --output {output} --mlst MLST_Nmen --amrfinder AMRfinder_Nmen --vfdb ABRicate_VFDB_Nmen {params.samples} 2>&1>{log}
    """

# Convert the summary file (tsv format) to an Excel file and save this version in the output folder
rule typing_summary_to_xlsx_Nmen:
  input:
    "backup/{timestamp}/Nmen_typing_summary.tsv"
  output:
    "output/{timestamp}/Nmen_typing_summary.xlsx"
  conda:
    "envs/python.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/Nmen_typing_summary_to_xlsx.log"
  threads: 16
  shell:
    """
    python workflow/scripts/summary_to_xlsx.py -d '\t' {input} {output} 2>&1>{log}
    """
