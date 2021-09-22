configfile: "workflow/config/config.yaml"

IDS, = glob_wildcards("input/{sample}_L001_R1_001.fastq.gz")

rule all:
  input:
    # "summary.xlsx",
    # expand("output/assemblies/{sample}.fasta", sample=IDS),
    "output/qc_reports",
    expand("backup/{sample}.tar.gz", sample=IDS)

rule fastqc_pre_R1:
  input:
    "input/{sample}_L001_R1_001.fastq.gz"
  output:
    html = "tmp_data/fastqc_pre_out/{sample}_R1.html",
    zip = "tmp_data/fastqc_pre_out/{sample}_R1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "logs/fastqc_pre/{sample}_R1.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

rule fastqc_pre_R2:
  input:
    "input/{sample}_L001_R2_001.fastq.gz"
  output:
    html = "tmp_data/fastqc_pre_out/{sample}_R2.html",
    zip = "tmp_data/fastqc_pre_out/{sample}_R2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "logs/fastqc_pre/{sample}_R2.log"
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

rule fastqc_post_R1:
  input:
    "tmp_data/trimmed/{sample}_L001_R1_001_corrected.fastq.gz",
  output:
    html = "tmp_data/fastqc_post_out/{sample}_R1.html",
    zip = "tmp_data/fastqc_post_out/{sample}_R1_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "logs/fastqc_post/{sample}_R1.log"
  threads: 2
  wrapper:
    "0.78.0/bio/fastqc"

rule fastqc_post_R2:
  input:
    "tmp_data/trimmed/{sample}_L001_R2_001_corrected.fastq.gz",
  output:
    html = "tmp_data/fastqc_post_out/{sample}_R2.html",
    zip = "tmp_data/fastqc_post_out/{sample}_R2_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file
  params: "--quiet"
  log:
    "logs/fastqc_post/{sample}_R2.log"
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

rule multiqc_fastqc:
  input:
    expand("tmp_data/fastqc_pre_out/{sample}_{read}_fastqc.zip", sample=IDS, read = ['R1', 'R2']),
    expand("tmp_data/fastqc_post_out/{sample}_{read}_fastqc.zip", sample=IDS, read = ['R1', 'R2']),
  output:
    "tmp_data/fastqc_report.html"
  log:
    "logs/multiqc_fastqc.log"
  wrapper:
    "0.78.0/bio/multiqc"

rule multiqc_kraken:
  input:
    expand("tmp_data/kraken_out/{sample}_kraken2_report.txt", sample=IDS),
  output:
    "tmp_data/kraken_report.html"
  log:
    "logs/multiqc_kraken.log"
  wrapper:
    "0.78.0/bio/multiqc"

rule multiqc_quast:
  input:
    expand("tmp_data/quast_out/{sample}", sample=IDS)
  output:
    "tmp_data/quast_report.html"
  log:
    "logs/multiqc_quast.log"
  wrapper:
    "0.78.0/bio/multiqc"

rule copy_data:
  input:
    "tmp_data/fastqc_report.html",
    "tmp_data/kraken_report.html",
    "tmp_data/quast_report.html"
  output:
    directory("output/qc_reports")
  threads: 1
  shell:
    """
    mkdir -p {output}
    cp {input} {output}
    """

rule backup_data:
  input:
    genome = "output/genomes/{sample}.fasta",
    kraken = "tmp_data/kraken_out/{sample}_kraken2_report.txt",
    quast = "tmp_data/quast_out/{sample}",
    fastqc_pre_R1 = "tmp_data/fastqc_pre_out/{sample}_R1.html",
    fastqc_pre_R2 = "tmp_data/fastqc_pre_out/{sample}_R2.html",
    fastqc_post_R1 = "tmp_data/fastqc_post_out/{sample}_R1.html",
    fastqc_post_R2 = "tmp_data/fastqc_post_out/{sample}_R2.html",
  output:
    "backup/{sample}.tar.gz"    
  params:
    sample = "{sample}"
  shell:
    """
    mkdir -p tmp_data/.backup/{params.sample}/quast
    mkdir -p tmp_data/.backup/{params.sample}/fastqc_pre
    mkdir -p tmp_data/.backup/{params.sample}/fastqc_post
    mkdir -p tmp_data/.backup/{params.sample}/kraken
    cp {input.fastqc_pre_R1} tmp_data/.backup/{params.sample}/fastqc_pre
    cp {input.fastqc_post_R1} tmp_data/.backup/{params.sample}/fastqc_post
    cp {input.fastqc_pre_R2} tmp_data/.backup/{params.sample}/fastqc_pre
    cp {input.fastqc_post_R2} tmp_data/.backup/{params.sample}/fastqc_post
    cp -r {input.quast} tmp_data/.backup/{params.sample}/quast
    cp {input.kraken} tmp_data/.backup/{params.sample}/kraken
    cp {input.genome} tmp_data/.backup/{params.sample}
    tar zcvf backup/{params.sample}.tar tmp_data/.backup/{params.sample}    
    pigz -9 backup/{params.sample}.tar
    """
