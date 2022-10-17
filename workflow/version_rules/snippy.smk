configfile: "workflow/config/config.yaml"

rule snippy_version:
  input:
    ref = "tmp_data/{timestamp}/references/MGAS5005.gbk"
  output:
    temp("tmp_data/{timestamp}/versions/snippy.txt")
  conda:
    "../envs/snippy.yaml"
  shell:
    """
    snippy --version | sed 's/^/snippy,/' > {output}
    snippy --check 2>&1 | grep 'Checking version' | cut -d " " -f 4,11-13 | sed 's/have //' | sed 's/ /,/' >> {output}
    echo {input.ref} | sed 's/^/snippy,reference_genome_name: /' >> {output}
    grep -E '^VERSION' {input.ref} | tr -s ' ' | cut -d ' ' -f 2 | sed 's/^/snippy,reference_genome_version: /' >> {output}
    md5sum {input.ref} | sed 's/^/snippy,reference_genome_md5sum: /' >> {output}
    stat {input.ref} | grep Modify | cut -f 2-4 -d " " | sed 's/^/snippy,reference_genome_date: /' >> {output}
    """
