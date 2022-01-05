# vim: set filetype=sh :

rule multiqc_quality_control_concatenatedData:
    """
    Quality control of concatenated fastq files (MULTIQC).
    """
    input:
        expand("%s/QC/CONCATENATED/{samples}_R1_fastqc.zip" % (config["project-folder"]), samples=samples)
    output:
        directory("%s/QC/CONCATENATED/multiqc/" % (config["project-folder"]))
    log:
        "%s/logs/MULTIQC/multiqc_concatenated.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/MULTIQC/multiqc_concatenated.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        multiqc -f -o {output} {input} 2> {log};
    """