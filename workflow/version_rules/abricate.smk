configfile: "workflow/config/config.yaml"

rule ABRicate_version:
  output:
    "tmp_data/{timestamp}/versions/ABRicate.txt"
  conda:
    "../envs/abricate.yaml"
  params:
    db_fimH = config["abricate"]["fimH_db"],
    db_vfdb = config["abricate"]["VFDB_db"], 
    db_datadir = config["abricate"]["fimH_datadir"],
  shell:
    """
    abricate --version | sed 's/^/ABRicate,/' > {output}
    blastn -version | head -n 1 | sed 's/^/blastn,/' >> {output}
    echo {params.db_fimH} | sed 's/^/ABRicate_fimH,db_name: /' >> {output}
    md5sum {params.db_datadir}/{params.db_fimH}/sequences | sed 's/^/ABRicate_fimH,db_md5sum: /' >> {output}
    stat {params.db_datadir}/{params.db_fimH}/sequences | grep Modify | cut -f 2-4 -d " " | sed 's/^/ABRicate_fimH,db_date: /' >> {output}
    ABRICATE_PATH=$(command -v abricate)
    echo {params.db_vfdb} | sed 's/^/ABRicate_VFDB,db_name: /' >> {output}
    md5sum ${{ABRICATE_PATH%/abricate}}/../db/{params.db_vfdb}/sequences | sed 's/^/ABRicate_VFDB,db_md5sum: /' >> {output}
    stat ${{ABRICATE_PATH%/abricate}}/../db/{params.db_vfdb}/sequences | grep Modify | cut -f 2-4 -d " " | sed 's/^/ABRicate_VFDB,db_date: /' >> {output}
    """
