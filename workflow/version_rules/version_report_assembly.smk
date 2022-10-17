configfile: "workflow/config/config.yaml"

rule version_report_assembly:
  input:
    general = "tmp_data/{timestamp}/versions/assembly_general.txt",
    coverage = "tmp_data/{timestamp}/versions/coverage.txt",
    kraken2 = "tmp_data/{timestamp}/versions/kraken.txt",
    mlst = "tmp_data/{timestamp}/versions/mlst.txt",
    quast = "tmp_data/{timestamp}/versions/quast.txt",
    shovill = "tmp_data/{timestamp}/shovill.txt",
  output:
    "backup/{timestamp}/version_report_assembly.csv"
  shell:
    """
    echo "rule,tool,version_data" > {output}
    sed 's/^/general,/' {input.general} >> {output}
    sed 's/^/coverage,/' {input.coverage} >> {output}
    sed 's/^/kraken2,/' {input.kraken2} >> {output}
    sed 's/^/mlst,/' {input.mlst} >> {output}
    sed 's/^/quast,/' {input.quast} >> {output}
    sed 's/^/shovill,/' {input.shovill} >> {output}
    """
