# vim: set filetype=sh :

rule fastqc_quality_control_ars_trimmed:
    """
    Quality control of fastq files (FASTQC).
    """
    input:
        unmapped="%s/FASTQ/ars_trimmed/unmapped/{samples}_ars_trimmed_unmapped.fastq" % (config["project-folder"]),
        mapped="%s/FASTQ/ars_trimmed/mapped/{samples}_ars_trimmed_mapped.fastq" % (config["project-folder"])
    output:
        unmapped="%s/QC/ars_trimmed/{samples}_ars_trimmed_unmapped_fastqc.zip" % (config["project-folder"]),
        mapped="%s/QC/ars_trimmed/{samples}_ars_trimmed_mapped_fastqc.zip" % (config["project-folder"])
    log:
        "%s/logs/FASTQC/fastqc_ars_trimmed.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/FASTQC/fastqc_ars_trimmed.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    params:
        outfolder="%s/QC/ars_trimmed/" % (config["project-folder"])
    shell:"""
        mkdir -p {params.outfolder};
        fastqc -t {threads} -o {params.outfolder} --extract {input.unmapped} 2> {log};
        fastqc -t {threads} -o {params.outfolder} --extract {input.mapped} 2> {log};
    """