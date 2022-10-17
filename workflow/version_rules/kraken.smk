configfile: "workflow/config/config.yaml"

rule kraken2_version:
  output:
    temp("tmp_data/{timestamp}/versions/kraken.txt")
  conda:
    "../envs/kraken.yaml"
  params:
    db = config["kraken"]["db"]
  shell:
    """
    kraken2 --version | head -n 1 | sed 's/^/kraken2,/' > {output}
    echo {params.db} | sed 's/^/kraken2,db_name: /' >> {output}
    md5sum {params.db}/hash.k2d | sed 's/^/kraken2,db_md5sum: /' >> {output}
    stat {params.db}/hash.k2d | grep Modify | cut -f 2-4 -d " " | sed 's/^/kraken2,db_date: /' >> {output}
    """
