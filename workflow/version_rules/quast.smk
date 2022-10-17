configfile: "workflow/config/config.yaml"

rule quast_version:
  output:
    temp("tmp_data/{timestamp}/versions/quast.txt")
  conda:
    "../envs/quast.yaml"
  shell:
    """
    quast --version | sed 's/^/quast,/' > {output}
    """
