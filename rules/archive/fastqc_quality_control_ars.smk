# vim: set filetype=sh :

rule fastqc_quality_control_ars:
    """
    Quality control of fastq files (FASTQC).
    """
    input:
        unmapped="%s/FASTQ/ars/unmapped/{samples}_ars_unmapped.fastq" % (config["project-folder"]),
        mapped="%s/FASTQ/ars/mapped/{samples}_ars_mapped.fastq" % (config["project-folder"])
    output:
        unmapped="%s/QC/ars/{samples}_ars_unmapped_fastqc.zip" % (config["project-folder"]),
        mapped="%s/QC/ars/{samples}_ars_mapped_fastqc.zip" % (config["project-folder"])
    log:
        "%s/logs/FASTQC/fastqc_ars.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/FASTQC/fastqc_ars.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    params:
        outfolder="%s/QC/ars/" % (config["project-folder"])
    shell:"""
        mkdir -p {params.outfolder};
        fastqc -t {threads} -o {params.outfolder} --extract {input.unmapped} 2> {log};
        fastqc -t {threads} -o {params.outfolder} --extract {input.mapped} 2> {log};
    """