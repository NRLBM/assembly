import sys

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

if len(SPYO_SAMPLES) > 0:
  all_output.append(expand("output/{timestamp}/Spyo_typing_summary.xlsx", timestamp=config["timestamp"]))
  all_output.append(expand("backup/{timestamp}/Spyo_typing_summary.tsv", timestamp=config["timestamp"]))
  all_output.append(expand("backup/{timestamp}/Spyo_typing/{sample}_typing.tar.gz", sample=SPYO_SAMPLES, timestamp=config["timestamp"]))

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
    amrfinder --threads 4 --nucleotide {input} --organism {params.organism} --output {output} 2>&1>{log}    
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
    expand("tmp_data/{timestamp}/ABRicate_VFDB_Ecoli/{sample}.tsv", sample=ECOLI_SAMPLES, timestamp=config["timestamp"]),
    qc_report = expand("backup/{timestamp}/qc_report_ecoli.tsv", timestamp=config["timestamp"]),
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
    python workflow/scripts/typing_summary.py --species Ecoli --timestamp {params.timestamp} --output {output} --qc {input.qc_report} --mlst MLST_Ecoli --amrfinder AMRfinder_Ecoli --vfdb ABRicate_VFDB_Ecoli {params.samples} 2>&1>{log}
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
    python workflow/scripts/bexsero_coverage.py --input {input.json} --output {output.tsv} 2>&1>{log}
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
    expand("tmp_data/{timestamp}/ABRicate_VFDB_Nmen/{sample}.tsv", sample=NMEN_SAMPLES, timestamp=config["timestamp"]),
    qc_report = expand("backup/{timestamp}/qc_report_nmen.tsv", timestamp=config["timestamp"]),
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
    python workflow/scripts/typing_summary.py --species Nmen --timestamp {params.timestamp} --qc {input.qc_report} --output {output} --mlst MLST_Nmen --amrfinder AMRfinder_Nmen --vfdb ABRicate_VFDB_Nmen {params.samples} 2>&1>{log}
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

# Streptococcus pyogenes typing
rule MLST_spyo:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta"
  output:
    json = temp("tmp_data/{timestamp}/MLST_Spyo/{sample}.json")
  threads: 8
  params:
    scheme = config['pubmlst']['spyo_mlst_scheme']
  conda:
    "envs/python.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/MLST_Spyo_PubMLST_API/{sample}.log"
  shell:
    """
    python workflow/scripts/type_pubmlst_api.py --genome {input.genome} --api-url {params.scheme} 1> {output.json} 2>{log}
    """

rule convert_MLST_spyo:
  input:
    "tmp_data/{timestamp}/MLST_Spyo/{sample}.json"
  output:
    temp("tmp_data/{timestamp}/MLST_Spyo/{sample}.tsv")
  threads: 1
  conda:
    "envs/python.yaml"
  params:
    loci = config['pubmlst']['spyo_mlst_loci']
  log:
    "slurm/snakemake_logs/{timestamp}/MLST_Spyo_convert/{sample}.log"
  shell:
    """
    python workflow/scripts/convert_MLST.py --input {input} --output {output} --loci {params.loci}
    """

rule emmtyper_spyo:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta"
  output:
    tsv = temp("tmp_data/{timestamp}/emmtyper_Spyo/{sample}.tsv")
  threads: 8
  conda:
    "envs/emmtyper.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/emmtyper_Spyo/{sample}.log"
  shell:
    """
    emmtyper --output-format verbose {input.genome} > {output.tsv}
    """

rule gunzip_ref:
  input:
    "workflow/references/MGAS5005.gbk.gz"
  output:
    temp("tmp_data/{timestamp}/references/MGAS5005.gbk")
  threads: 16
  shell:
    """
    gunzip -c {input} > {output}
    """

rule snippy_M1UK_Spyo:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta",
    ref = "tmp_data/{timestamp}/references/MGAS5005.gbk"
  output:
    temp(directory("tmp_data/{timestamp}/snippy_M1UK_Spyo/{sample}"))
  log:
    "slurm/snakemake_logs/{timestamp}/snippy_M1UK_Spyo/{sample}.log"
  params:
    general = config["snippy"]["general"]
  threads: 16
  conda: "envs/snippy.yaml"
  shell:
    """
    snippy {params.general} --cpus {threads} --outdir {output} --ref {input.ref} --contigs {input.genome} 2>&1>{log}
    """

rule compare_SNPs:
  input:
    snippy = "tmp_data/{timestamp}/snippy_M1UK_Spyo/{sample}",
    ref = "workflow/references/M1UK.vcf"
  output:
    "tmp_data/{timestamp}/M1UK_SNPs_Spyo/{sample}.tsv"
  conda:
    "envs/compare_SNPs.yaml"
  threads: 1
  log:
    "slurm/snakemake_logs/{timestamp}/compare_SNPs/{sample}.log"
  shell:
    """
    bgzip --keep --force {input.snippy}/snps.vcf 2>&1>{log}
    bcftools index {input.snippy}/snps.vcf.gz 2>&1>>{log}
    bcftools view -O v -R {input.ref} {input.snippy}/snps.vcf.gz 2>>{log} | vt decompose_blocksub -o - - 2>>{log} | python workflow/scripts/compare_SNPs.py --ref {input.ref} --extended > {output} 2>>{log}
    """

# Find all files to backup and create a tarball with a timestamp
rule backup_data_spyo:
  input:
    M1UK_tsv = "tmp_data/{timestamp}/M1UK_SNPs_Spyo/{sample}.tsv",
    emmtyper_tsv = "tmp_data/{timestamp}/emmtyper_Spyo/{sample}.tsv",
    MLST_tsv = "tmp_data/{timestamp}/MLST_Spyo/{sample}.tsv",
    MLST_json = "tmp_data/{timestamp}/MLST_Spyo/{sample}.json",
  output:
    "backup/{timestamp}/Spyo_typing/{sample}_typing.tar.gz"
  params:
    sample = "{sample}",
    timestamp = config["timestamp"]
  threads: 16
  shell:
    """
    mkdir -p tmp_data/.backup/{params.sample}/MLST_PubMLST
    mkdir -p tmp_data/.backup/{params.sample}/emmtyper
    mkdir -p tmp_data/.backup/{params.sample}/M1UK_SNPs
    cp {input.emmtyper_tsv} tmp_data/.backup/{params.sample}/emmtyper
    cp {input.MLST_tsv} tmp_data/.backup/{params.sample}/MLST_PubMLST
    cp {input.MLST_json} tmp_data/.backup/{params.sample}/MLST_PubMLST
    cp {input.M1UK_tsv} tmp_data/.backup/{params.sample}/M1UK_SNPs
    tar zcvf {output} -C tmp_data/.backup {params.sample}
    rm -rf tmp_data/.backup
    """

rule typing_summary_csv_Spyo:
  input:
    expand("tmp_data/{timestamp}/M1UK_SNPs_Spyo/{sample}.tsv", sample=SPYO_SAMPLES, timestamp=config["timestamp"]),
    expand("tmp_data/{timestamp}/emmtyper_Spyo/{sample}.tsv", sample=SPYO_SAMPLES, timestamp=config["timestamp"]),
    expand("tmp_data/{timestamp}/MLST_Spyo/{sample}.tsv", sample=SPYO_SAMPLES, timestamp=config["timestamp"]),
    qc_report = expand("backup/{timestamp}/qc_report_spyo.tsv", timestamp=config["timestamp"]),
  output:
    "backup/{timestamp}/Spyo_typing_summary.tsv"
  conda:
    "envs/python.yaml"
  params:
    timestamp = config["timestamp"],
    samples = SPYO_SAMPLES
  threads: 16
  log:
    "slurm/snakemake_logs/{timestamp}/Spyo_typing_summary.log"
  shell:
    """
    python workflow/scripts/typing_summary.py --species Spyo --timestamp {params.timestamp} --qc {input.qc_report} --output {output} --M1UK M1UK_SNPs_Spyo --mlst MLST_Spyo --emmtyper emmtyper_Spyo {params.samples} 2>&1>{log}
    """

# Convert the summary file (tsv format) to an Excel file and save this version in the output folder
rule typing_summary_to_xlsx_Spyo:
  input:
    "backup/{timestamp}/Spyo_typing_summary.tsv"
  output:
    "output/{timestamp}/Spyo_typing_summary.xlsx"
  conda:
    "envs/python.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/Spyo_typing_summary_to_xlsx.log"
  threads: 16
  shell:
    """
    python workflow/scripts/summary_to_xlsx.py -d '\t' {input} {output} 2>&1>{log}
    """
