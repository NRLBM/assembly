configfile: "workflow/config/config.yaml"

rule AMRFinder_version:
  output: 
    temp("tmp_data/{timestamp}/versions/AMRfinder.txt")
  conda:
    "../envs/amrfinder.yaml"
  shell:
    """
    amrfinder --version | sed 's/^/AMRfinder,/' > {output}
    AMRFINDER_PATH=$(command -v amrfinder)
    cat ${{AMRFINDER_PATH%/amrfinder}}/../share/amrfinderplus/data/latest/version.txt | sed 's/^/AMRfinder,db_version: /' >> {output}
    """
