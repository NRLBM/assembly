configfile: "workflow/config/config.yaml"

rule mlst_version:
  output:
    temp("tmp_data/{timestamp}/versions/mlst.txt")
  conda:
    "../envs/mlst.yaml"
  shell:
    """
    mlst --version | sed 's/^/mlst,/' > {output}
    MLST_PATH=$(command -v mlst)
    stat ${{MLST_PATH%/mlst}}/../db/blast/mlst.fa | grep Modify | cut -f 2-4 -d " " | sed 's/^/mlst,db_date: /' >> {output}
    """
