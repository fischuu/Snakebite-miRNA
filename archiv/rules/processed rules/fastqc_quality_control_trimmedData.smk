# vim: set filetype=sh :

rule fastqc_quality_control_trimmedData:
    """
    Quality control of fastq files (FASTQC).
    """
    input:
        "%s/FASTQ/TRIMMED/{raw_samples}_trimmed.fastq.gz" % (config["project-folder"])
    output:
        "%s/QC/TRIMMED/{raw_samples}_trimmed_fastqc.zip" % (config["project-folder"])
    log:
        "%s/logs/FASTQC/fastqc_trimmed.{raw_samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/FASTQC/fastqc_trimmed.{raw_samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    params:
        outfolder="%s/QC/TRIMMED/" % (config["project-folder"])
    shell:"""
        mkdir -p {params.outfolder};
        fastqc -t {threads} -o {params.outfolder} --extract {input} 2> {log};
    """