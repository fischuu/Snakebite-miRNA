# vim: set filetype=sh :

rule multiqc_quality_control_mergedData_ars_trimmed:
    """
    Quality control of fastq files (MULTIQC).
    """
    input:
        unmapped=expand("%s/QC/ars_trimmed/{samples}_ars_trimmed_unmapped_fastqc.zip" % (config["project-folder"]), samples=samples),
        mapped=expand("%s/QC/ars_trimmed/{samples}_ars_trimmed_mapped_fastqc.zip" % (config["project-folder"]), samples=samples)
    output:
        unmapped=directory("%s/QC/ars_trimmed/multiqc_unmapped/" % (config["project-folder"])),
        mapped=directory("%s/QC/ars_trimmed/multiqc_mapped/" % (config["project-folder"]))
    log:
        "%s/logs/MULTIQC/multiqc_ars_trimmed.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/MULTIQC/multiqc_ars_trimmed.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        multiqc -f -o {output.unmapped} {input.unmapped} 2> {log};
        multiqc -f -o {output.mapped} {input.mapped} 2> {log};
    """