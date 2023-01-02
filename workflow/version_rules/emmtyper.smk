configfile: "workflow/config/config.yaml"

rule emmtyper_version:
  input:
    dummy = "tmp_data/{timestamp}/emmtyper_db_updated.txt",
  output:
    temp("tmp_data/{timestamp}/versions/emmtyper.txt")
  params:
    db = config["emmtyper"]["db"]
  conda:
    "../envs/emmtyper.yaml"
  shell:
    """
    emmtyper --version | sed 's/^/emmtyper,/' > {output}
    echo {params.db} | sed 's/^/emmtyper,db_name: /' >> {output}
    md5sum {params.db} | sed 's/^/emmtyper,db_md5sum: /' >> {output}
    stat {params.db} | grep Modify | cut -f 2-4 -d " " | sed 's/^/emmtyper,db_date: /' >> {output}
    """
