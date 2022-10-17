configfile: "workflow/config/config.yaml"

rule combine_mash_nmen_version:
  output:
    temp("tmp_data/{timestamp}/versions/combine_mash_nmen.txt")
  conda:
    "../envs/combine_mash_Nmen.yaml"
  shell:
    """
    python workflow/scripts/combine_mash.py --version | sed 's/^/combine_mash.py,/' > {output}
    python -c "import pandas; print(pandas.__version__)" | sed 's/^/pandas,/' >> {output}
    python -c "import networkx; print(networkx.__version__)" | sed 's/^/networkx,/' >> {output}
    """

