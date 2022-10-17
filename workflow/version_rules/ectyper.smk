configfile: "workflow/config/config.yaml"

rule ECtyper_version:
  output:
    temp("tmp_data/{timestamp}/versions/ectyper.txt")
  conda:
    "../envs/ectyper.yaml"
  shell:
    """
    ectyper --version | sed 's/^/ectyper,/' > {output}
    """
