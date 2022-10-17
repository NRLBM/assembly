configfile: "workflow/config/config.yaml"

rule meningotype_version:
  output:
    "tmp_data/{timestamp}/versions/meningotype.txt"
  conda:
    "../envs/meningotype.yaml"
  shell:
    """
    meningotype --version | sed 's/^/meningotype,/' > {output}
    MENINGOTYPE_PATH=$(command -v meningotype)
    stat ${{MENINGOTYPE_PATH%/meningotype}}/../lib/python3.10/site-packages/meningotype/db/PorA_VR1.fas | grep Modify | cut -f 2-4 -d " " | sed 's/^/meningotype,db_PorA_VR1_date: /' >> {output}
    stat ${{MENINGOTYPE_PATH%/meningotype}}/../lib/python3.10/site-packages/meningotype/db/PorA_VR2.fas | grep Modify | cut -f 2-4 -d " " | sed 's/^/meningotype,db_PorA_VR2_date: /' >> {output}
    stat ${{MENINGOTYPE_PATH%/meningotype}}/../lib/python3.10/site-packages/meningotype/db/FetA_VR.fas | grep Modify | cut -f 2-4 -d " " | sed 's/^/meningotype,db_FetA_VR_date: /' >> {output}
    """
