configfile: "workflow/config/config.yaml"

rule assembly_general_version:
  output:
    temp("tmp_data/{timestamp}/versions/assembly_general.txt"),
  params:
    wrapper_trimmomatic_pe = config["wrapper"]["trimmomatic_pe"],
    wrapper_fastqc = config["wrapper"]["fastqc"],
    wrapper_multiqc = config["wrapper"]["multiqc"]
  conda:
    "../envs/python.yaml"
  shell:
    """
    python workflow/scripts/move_to_scratch.py --version | sed 's/^/move_to_scratch.py,/' > {output}
    python workflow/scripts/summary.py --version | sed 's/^/summary.py,/' >> {output}
    python workflow/scripts/summary_to_xlsx.py --version | sed 's/^/summary_to_xlsx.py,/' >> {output}
    python -c "import pandas; print(pandas.__version__)" | sed 's/^/pandas,/' >> {output}
    python -c "import openpyxl; print(openpyxl.__version__)" | sed 's/^/openpyxl,/' >> {output}
    echo {params.wrapper_trimmomatic_pe} | sed 's/^/trimmomatic_pe,snakemake_wrapper: /' >> {output}
    echo {params.wrapper_fastqc} | sed 's/^/fastqc,snakemake_wrapper: /' >> {output}
    echo {params.wrapper_multiqc}  | sed 's/^/multiqc,snakemake_wrapper: /' >> {output}
    """
