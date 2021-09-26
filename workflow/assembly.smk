configfile: "workflow/config/config.yaml"

IDS, = glob_wildcards("input/{sample}_L001_R1_001.fastq.gz")

rule all:
  input:
    expand("output/{timestamp}/qc_reports", timestamp=config["timestamp"]),
    expand("backup/{timestamp}/{sample}.tar.gz", sample=IDS, timestamp=config["timestamp"]),
    expand("output/{timestamp}/summary.xlsx", timestamp=config["timestamp"])

rule fastqc_pre_R1:
  input:
    "input/{sample}_L001_R1_001.fastq.gz"
  output:
    html = temp("tmp_data/{timestamp}/fastqc_pre_out/{sample}_R1.html"),
    zip = temp("tmp_data/{timestamp}/fastqc_pre_out/{sample}_R1_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "slurm/snakemake_logs/{timestamp}/fastqc_pre/{sample}_R1.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

rule fastqc_pre_R2:
  input:
    "input/{sample}_L001_R2_001.fastq.gz"
  output:
    html = temp("tmp_data/{timestamp}/fastqc_pre_out/{sample}_R2.html"),
    zip = temp("tmp_data/{timestamp}/fastqc_pre_out/{sample}_R2_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "slurm/snakemake_logs/{timestamp}/fastqc_pre/{sample}_R2.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

rule trimmomatic_pe:
  input:
    r1 = "input/{sample}_L001_R1_001.fastq.gz",
    r2 = "input/{sample}_L001_R2_001.fastq.gz"
  output:
    r1 = temp("tmp_data/{timestamp}/trimmed/{sample}_L001_R1_001_corrected.fastq.gz"),
    r2 = temp("tmp_data/{timestamp}/trimmed/{sample}_L001_R2_001_corrected.fastq.gz"),
    r1_unpaired = temp("tmp_data/{timestamp}/trimmed/{sample}.1.unpaired.fastq.gz"),
    r2_unpaired = temp("tmp_data/{timestamp}/trimmed/{sample}.2.unpaired.fastq.gz")
  log:
    "slurm/snakemake_logs/{timestamp}/trimmomatic/{sample}.log"
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

rule fastqc_post_R1:
  input:
    "tmp_data/{timestamp}/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
  output:
    html = temp("tmp_data/{timestamp}/fastqc_post_out/{sample}_R1.html"),
    zip = temp("tmp_data/{timestamp}/fastqc_post_out/{sample}_R1_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "slurm/snakemake_logs/{timestamp}/fastqc_post/{sample}_R1.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

rule fastqc_post_R2:
  input:
    "tmp_data/{timestamp}/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
  output:
    html = temp("tmp_data/{timestamp}/fastqc_post_out/{sample}_R2.html"),
    zip = temp("tmp_data/{timestamp}/fastqc_post_out/{sample}_R2_fastqc.zip") # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "slurm/snakemake_logs/{timestamp}/fastqc_post/{sample}_R2.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

rule shovill:
  input:
    fw = "tmp_data/{timestamp}/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    rv = "tmp_data/{timestamp}/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
  output:
    assembly = "output/{timestamp}/genomes/{sample}.fasta",
    shovill = temp(directory("tmp_data/{timestamp}/shovill_out/{timestamp}_{sample}")),
  conda:
    "envs/shovill.yaml"
  params:
    minlen = config["shovill"]["minlen"],
    ram = config["shovill"]["ram"],
    depth = config["shovill"]["depth"],
    assembler = config["shovill"]["assembler"],
    tmpdir = config["shovill"]["tmpdir"]
  log:
    "slurm/snakemake_logs/{timestamp}/shovill/{sample}.log"
  threads: 16
  shell:
    """
    GENOME_SIZE=$(bash workflow/scripts/estimated_genome_size.sh {input.fw})
    shovill --assembler {params.assembler} --outdir {output.shovill} --tmp {params.tmpdir} --depth {params.depth} --gsize ${{GENOME_SIZE}} --cpus {threads} --ram {params.ram} --minlen {params.minlen} --R1 {input.fw} --R2 {input.rv} 2>&1>{log}
    cp {output.shovill}/contigs.fa {output.assembly}
    """

rule coverage:
  input:
    fw = "tmp_data/{timestamp}/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    rv = "tmp_data/{timestamp}/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
    assembly = "output/{timestamp}/genomes/{sample}.fasta",
  output:
    temp("tmp_data/{timestamp}/coverage_out/{sample}.txt")
  conda:
    "envs/coverage.yaml"
  params:
    minimap_x = "sr"
  log:
    "slurm/snakemake_logs/{timestamp}/coverage/{sample}.log"
  threads: 16
  shell:
    """
    minimap2 -a -x {params.minimap_x} -t {threads} {input.assembly} {input.fw} {input.rv} | samtools sort -l 0 --threads {threads} | bedtools genomecov -d -ibam stdin | awk '{{t += $3}} END {{print t/NR}}' 1>{output} 2>{log}
    """

rule quast:
  input:
    "output/{timestamp}/genomes/{sample}.fasta"
  output:
    temp(directory("tmp_data/{timestamp}/quast_out/{sample}"))
  conda:
    "envs/quast.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/quast/{sample}.log"
  threads: 8
  shell:
    """
    quast --threads {threads} -o {output} {input} 2>&1>{log}
    """

rule kraken2:
  input:
    fw = "tmp_data/{timestamp}/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
    rv = "tmp_data/{timestamp}/trimmed/{sample}_L001_R2_001_corrected.fastq.gz"
  output:
    report = temp("tmp_data/{timestamp}/kraken_out/{sample}.txt")
  conda:
    "envs/kraken.yaml"
  params:
    general = config["kraken"]["general"],
    db = config["kraken"]["db"]
  log:
    "slurm/snakemake_logs/{timestamp}/kraken2/{sample}.log"
  threads: 8
  shell:
    """
    kraken2 --db {params.db} {params.general} --threads {threads} --report {output.report} {input.fw} {input.rv} 2>&1>{log}
    """

rule mlst:
  input:
    "output/{timestamp}/genomes/{sample}.fasta"
  output:
    temp("tmp_data/{timestamp}/mlst_out/{sample}.txt")
  log:
    "slurm/snakemake_logs/{timestamp}/mlst/{sample}/log"
  conda:
    "envs/mlst.yaml"
  threads: 1
  shell:
    """
    mlst {input} 1>{output} 2>{log}
    """

rule multiqc_fastqc:
  input:
    expand("tmp_data/{timestamp}/fastqc_pre_out/{sample}_{read}_fastqc.zip", sample=IDS, read = ['R1', 'R2'], timestamp=config["timestamp"]),
    expand("tmp_data/{timestamp}/fastqc_post_out/{sample}_{read}_fastqc.zip", sample=IDS, read = ['R1', 'R2'], timestamp=config["timestamp"]),
  output:
    temp("tmp_data/{timestamp}/fastqc_report.html")
  log:
    "slurm/snakemake_logs/{timestamp}/multiqc_fastqc.log"
  wrapper:
    "0.78.0/bio/multiqc"

rule multiqc_kraken:
  input:
    expand("tmp_data/{timestamp}/kraken_out/{sample}.txt", sample=IDS, timestamp=config["timestamp"]),
    expand("tmp_data/{timestamp}/quast_out/{sample}", sample=IDS, timestamp=config["timestamp"])
  output:
    temp("tmp_data/{timestamp}/kraken_quast_report.html")
  log:
    "slurm/snakemake_logs/{timestamp}/multiqc_kraken_quast.log"
  wrapper:
    "0.78.0/bio/multiqc"

rule copy_qc_reports:
  input:
    "tmp_data/{timestamp}/fastqc_report.html",
    "tmp_data/{timestamp}/kraken_quast_report.html",
  output:
    directory("output/{timestamp}/qc_reports")
  threads: 1
  shell:
    """
    mkdir -p {output}
    cp {input} {output}
    """

rule backup_data:
  input:
    genome = "output/{timestamp}/genomes/{sample}.fasta",
    kraken = "tmp_data/{timestamp}/kraken_out/{sample}.txt",
    quast = "tmp_data/{timestamp}/quast_out/{sample}",
    fastqc_pre_R1 = "tmp_data/{timestamp}/fastqc_pre_out/{sample}_R1.html",
    fastqc_pre_R2 = "tmp_data/{timestamp}/fastqc_pre_out/{sample}_R2.html",
    fastqc_post_R1 = "tmp_data/{timestamp}/fastqc_post_out/{sample}_R1.html",
    fastqc_post_R2 = "tmp_data/{timestamp}/fastqc_post_out/{sample}_R2.html",
    mlst = "tmp_data/{timestamp}/mlst_out/{sample}.txt",
    coverage = "tmp_data/{timestamp}/coverage_out/{sample}.txt",
  output:
    "backup/{timestamp}/{sample}.tar.gz"    
  params:
    sample = "{sample}",
    timestamp = config["timestamp"]
  shell:
    """
    mkdir -p tmp_data/.backup/{params.sample}/quast
    mkdir -p tmp_data/.backup/{params.sample}/fastqc_pre
    mkdir -p tmp_data/.backup/{params.sample}/fastqc_post
    mkdir -p tmp_data/.backup/{params.sample}/kraken
    mkdir -p tmp_data/.backup/{params.sample}/coverage
    mkdir -p tmp_data/.backup/{params.sample}/mlst
    cp {input.fastqc_pre_R1} tmp_data/.backup/{params.sample}/fastqc_pre
    cp {input.fastqc_post_R1} tmp_data/.backup/{params.sample}/fastqc_post
    cp {input.fastqc_pre_R2} tmp_data/.backup/{params.sample}/fastqc_pre
    cp {input.fastqc_post_R2} tmp_data/.backup/{params.sample}/fastqc_post
    cp -r {input.quast} tmp_data/.backup/{params.sample}/quast
    cp {input.kraken} tmp_data/.backup/{params.sample}/kraken
    cp {input.genome} tmp_data/.backup/{params.sample}
    cp {input.mlst} tmp_data/.backup/{params.sample}/mlst
    cp {input.coverage} tmp_data/.backup/{params.sample}/coverage
    tar zcvf backup/{params.timestamp}/{params.sample}.tar.gz -C tmp_data/.backup {params.sample}
    rm -rf tmp_data/.backup
    """

rule summary:
  input:
    kraken = expand("tmp_data/{timestamp}/kraken_out/{sample}.txt", sample=IDS, timestamp=config["timestamp"]),
    quast = expand("tmp_data/{timestamp}/quast_out/{sample}", sample=IDS, timestamp=config["timestamp"]),
    fastqc_pre = expand("tmp_data/{timestamp}/fastqc_pre_out/{sample}_{read}_fastqc.zip", sample=IDS, read = ['R1', 'R2'], timestamp=config["timestamp"]),
    fastqc_post = expand("tmp_data/{timestamp}/fastqc_post_out/{sample}_{read}_fastqc.zip", sample=IDS, read = ['R1', 'R2'], timestamp=config["timestamp"]),
    mlst = expand("tmp_data/{timestamp}/mlst_out/{sample}.txt", sample=IDS, timestamp=config["timestamp"]),
    coverage = expand("tmp_data/{timestamp}/coverage_out/{sample}.txt", sample=IDS, timestamp=config["timestamp"]),
  output:
    "backup/{timestamp}/summary.csv"
  conda:
    "envs/python.yaml"
  params:
    timestamp = config["timestamp"],
    samples = IDS
  log:
    "slurm/snakemake_logs/{timestamp}/summary.log"
  threads: 1
  shell:
    """
    python workflow/scripts/summary.py --timestamp {params.timestamp} {params.samples} 1> {output} 2>{log}
    """

rule summary_to_xlsx:
  input:
    "backup/{timestamp}/summary.csv"
  output:
    "output/{timestamp}/summary.xlsx"
  conda:
    "envs/python.yaml"
  log:
    "slurm/snakemake_logs/{timestamp}/summary_to_xlsx.log"
  threads: 1
  shell:
    """
    python workflow/scripts/summary_to_xlsx.py {input} {output} 2>&1>{log}
    """
