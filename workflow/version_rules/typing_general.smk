configfile: "workflow/config/config.yaml"

rule typing_general_version:
  output:
    "tmp_data/{timestamp}/versions/typing_general.txt"
  conda:
    "../envs/python.yaml"
  shell:
    """
    python workflow/scripts/bexsero_coverage.py --version | sed 's/^/bexsero_coverage.py,/' > {output}
    python workflow/scripts/convert_MLST.py --version | sed 's/^/convert_MLST.py,/' >> {output}
    python workflow/scripts/summary_to_xlsx.py --version | sed 's/^/summary_to_xlsx.py,/' >> {output}
    python workflow/scripts/type_pubmlst_api.py --version | sed 's/^/type_pubmlst_api.py,/' >> {output}
    python workflow/scripts/typing_summary.py --version | sed 's/^/typing_summary.py,/' >> {output}
    """
