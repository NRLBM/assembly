rule qc_version:
  output:
    temp("tmp_data/{timestamp}/versions/qc.txt"),
  shell:
    """
    python workflow/scripts/qc.py --version | sed 's/^/qc,/' > {output}
    """
