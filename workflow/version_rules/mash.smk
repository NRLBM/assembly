configfile: "workflow/config/config.yaml"

rule mash_version:
  output:
    temp("tmp_data/{timestamp}/versions/mash.txt")
  conda:
    "../envs/mash.yaml"
  shell:
    """
    mash --version | sed 's/^/mash,/' > {output}
    """

