configfile: "workflow/config/config.yaml"

IDS, = glob_wildcards("input/{sample}_L001_R1_001.fastq.gz")

rule all:
  input:
    # "summary.xlsx",
    # expand("output/assemblies/{sample}.fasta", sample=IDS),
    "output/qc_report.html"

rule fastqc_pre:
  input:
    fw = "input/{sample}_L001_R1_001.fastq.gz",
    rv = "input/{sample}_L001_R1_001.fastq.gz"
  output:
    html = "tmp_data/fastqc_pre_out/{sample}.html",
    zip = "tmp_data/fastqc_pre_out/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "logs/fastqc_pre/{sample}.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

rule trimmomatic_pe:
  input:
    r1 = "input/{sample}_L001_R1_001.fastq.gz",
    r2 = "input/{sample}_L001_R2_001.fastq.gz"
  output:
    r1 = "tmp_data/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    r2 = "tmp_data/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
    # reads where trimming entirely removed the mate
    r1_unpaired="tmp_data/trimmed/{sample}.1.unpaired.fastq.gz",
    r2_unpaired="tmp_data/trimmed/{sample}.2.unpaired.fastq.gz"
  log:
    "logs/trimmomatic/{sample}.log"
  params:
    # list of trimmers (see manual)
    trimmer=["ILLUMINACLIP:all_paired.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36"],
    # optional parameters
    extra="",
    compression_level="-9"
  threads:
    16
  wrapper:
    "0.78.0/bio/trimmomatic/pe"

rule fastqc_post:
  input:
    r1 = "tmp_data/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    r2 = "tmp_data/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
  output:
    html = "tmp_data/fastqc_post_out/{sample}.html",
    zip = "tmp_data/fastqc_post_out/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "logs/fastqc_post/{sample}.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

rule shovill:
  input:
    fw = "tmp_data/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    rv = "tmp_data/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
  output:
    assembly = "output/genomes/{sample}.fasta",
    shovill = directory("tmp_data/shovill_out/{sample}"),
  conda:
    "envs/shovill.yaml"
  params:
    minlen = config["shovill"]["minlen"],
    ram = config["shovill"]["ram"],
    depth = config["shovill"]["depth"],
    assembler = config["shovill"]["assembler"],
    tmpdir = config["shovill"]["tmpdir"]
  log:
    "logs/shovill/{sample}.log"
  threads: 16
  shell:
    """
    GENOME_SIZE=$(bash workflow/scripts/estimated_genome_size.sh {input.fw})
    shovill --assembler {params.assembler} --outdir {output.shovill} --tmp {params.tmpdir} --depth {params.depth} --gsize ${{GENOME_SIZE}} --cpus {threads} --ram {params.ram} --minlen {params.minlen} --R1 {input.fw} --R2 {input.rv} 2>&1>{log}
    cp {output.shovill}/contigs.fa {output.assembly}
    """

rule quast:
  input:
    "output/genomes/{sample}.fasta"
  output:
    directory("tmp_data/quast_out/{sample}")
  conda:
    "envs/quast.yaml"
  log:
    "logs/quast/{sample}.log"
  threads: 8
  shell:
    """
    quast --threads {threads} -o {output} {input} 2>&1>{log}
    """

rule kraken2:
  input:
    fw = "tmp_data/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    rv = "tmp_data/trimmed/{sample}_L001_R2_001_corrected.fastq.gz"
  output:
    report = "tmp_data/kraken_out/{sample}_kraken2_report.txt"
  conda:
    "envs/kraken.yaml"
  params:
    general = config["kraken"]["general"],
    db = config["kraken"]["db"]
  log:
    "logs/kraken2/{sample}.log"
  threads: 8
  shell:
    """
    kraken2 --db {params.db} {params.general} --threads {threads} --report {output.report} {input.fw} {input.rv} 2>&1>{log}
    """

rule multiqc:
  input:
    expand("tmp_data/fastqc_pre_out/{sample}_fastqc.zip", sample=IDS),
    expand("tmp_data/fastqc_post_out/{sample}_fastqc.zip", sample=IDS),
    expand("tmp_data/kraken_out/{sample}_kraken2_report.txt", sample=IDS),
    expand("tmp_data/quast_out/{sample}", sample=IDS)
  output:
    "output/qc_report.html"
  log:
    "logs/multiqc.log"
  wrapper:
    "0.78.0/bio/multiqc"
