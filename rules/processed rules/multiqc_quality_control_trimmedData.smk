# vim: set filetype=sh :

rule multiqc_quality_control_trimmedData:
    """
    Quality control of fastq files (MULTIQC).
    """
    input:
        expand("%s/QC/TRIMMED/{temp2}_trimmed_fastqc.zip" % (config["project-folder"]), temp2=rawsamples)
    output:
        directory("%s/QC/TRIMMED/multiqc/" % (config["project-folder"]))
    log:
        "%s/logs/MULTIQC/multiqc_trimmed.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/MULTIQC/multiqc_trimmed.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        multiqc -f -o {output} {input} 2> {log};
    """