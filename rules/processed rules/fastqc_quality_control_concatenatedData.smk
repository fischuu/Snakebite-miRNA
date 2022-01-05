# vim: set filetype=sh :

rule fastqc_quality_control_concatenatedData:
    """
    Quality control of concatenated fastq files (FASTQC).
    """
    input:
        "%s/FASTQ/CONCATENATED/{samples}_R1.fastq.gz" % (config["project-folder"])
    output:
        "%s/QC/CONCATENATED/{samples}_R1_fastqc.zip" % (config["project-folder"])
    log:
        "%s/logs/FASTQC/fastqc_concatenated.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/FASTQC/fastqc_concatenated.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    params:
        outfolder="%s/QC/CONCATENATED/" % (config["project-folder"])
    shell:"""
        mkdir -p {params.outfolder};
        fastqc -t {threads} -o {params.outfolder} --extract {input} 2> {log};
    """