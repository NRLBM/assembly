configfile: "workflow/config/config.yaml"

rule coverage_version:
  output:
    temp("tmp_data/{timestamp}/versions/coverage.txt")
  conda:
    "../envs/coverage.yaml"
  shell:
    """
    minimap2 --version | sed 's/^/minimap2,/' > {output}
    samtools --version | head -n 2 | sed 's/^/samtools,/' >> {output}
    bedtools --version | sed 's/^/bedtools,/' >> {output}
    """

