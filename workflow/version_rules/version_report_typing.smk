configfile: "workflow/config/config.yaml"

rule version_report_assembly:
  input:
    qc = "tmp_data/{timestamp}/versions/qc.txt",
    general = "tmp_data/{timestamp}/versions/typing_general.txt",
    abricate = "tmp_data/{timestamp}/versions/ABRicate.txt",
    combine_mash_nmen = "tmp_data/{timestamp}/versions/combine_mash_nmen.txt",
    compare_SNPs = "tmp_data/{timestamp}/versions/compare_SNPs.txt",
    ectyper = "tmp_data/{timestamp}/versions/ectyper.txt",
    emmtyper = "tmp_data/{timestamp}/versions/emmtyper.txt",
    mash = "tmp_data/{timestamp}/versions/mash.txt",
    amrfinder = "tmp_data/{timestamp}/versions/AMRfinder.txt",
    meningotype = "tmp_data/{timestamp}/versions/meningotype.txt",
    snippy = "tmp_data/{timestamp}/versions/snippy.txt",
  output:
    "backup/{timestamp}/version_report_typing.csv"
  shell:
    """
    echo "rule,tool,version_data" > {output}
    sed 's/^/qc,/' {input.qc} >> {output}
    sed 's/^/general,/' {input.general} >> {output}
    sed 's/^/ABRicate,/' {input.abricate} >> {output}
    sed 's/^/combine_mash_nmen,/' {input.combine_mash_nmen} >> {output}
    sed 's/^/compare_SNPs,/' {input.compare_SNPs} >> {output}
    sed 's/^/ectyper,/' {input.ectyper} >> {output}
    sed 's/^/emmtyper,/' {input.emmtyper} >> {output}
    sed 's/^/mash,/' {input.mash} >> {output}
    sed 's/^/amrfinder,/' {input.amrfinder} >> {output}
    sed 's/^/meningotype,/' {input.meningotype} >> {output}
    sed 's/^/snippy,/' {input.snippy} >> {output}
    """
