configfile: "workflow/config/config.yaml"

rule compare_SNPs_version:
  input:
    "workflow/references/M1UK.vcf"
  output:
    versions = temp("tmp_data/{timestamp}/versions/compare_SNPs.txt"),
    tmp_vt = temp("tmp_data/{timestamp}/versions/tmp_vt.txt")
  conda:
    "../envs/compare_SNPs.yaml"
  shell:
    """
    bgzip --help | head -n 2 | tail -n 1 | sed 's/^/bgzip,/' > {output.versions}
    echo {input} | sed 's/^/M1UK_VCF,name: /' >> {output.versions}
    md5sum {input} | sed 's/^/M1UK_VCF,md5sum: /' >> {output.versions}
    stat {input} | grep Modify | cut -f 2-4 -d " " | sed 's/^/M1UK_VCF,date: /' >> {output.versions}
    bcftools --version-only | sed 's/^/bcftools,/' >> {output.versions}
    # For some reason, redirecting vt version into STDOUT and piping it gives an error. Using intermediate file remedies this 
    vt decompose -? 2>{output.tmp_vt}
    head -n 1 {output.tmp_vt} | sed 's/^/vt,/' >> {output.versions}
    python workflow/scripts/compare_SNPs.py --version | sed 's/^/compare_SNPs.py,/' >> {output.versions}
    """
