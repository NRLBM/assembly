configfile: "workflow/config/config.yaml"

rule shovill_version:
  output:
    temp("tmp_data/{timestamp}/shovill.txt"),
  conda:
    "../envs/shovill.yaml"
  shell:
    """
    shovill --version | sed 's/^/shovill,/' > {output}
    shovill --check 2>&1 | cut -f 3,7-99 -d " " | sed 's/ /,/' >> {output}
    mash --version | sed 's/^/mash,/' >> {output}
    """
