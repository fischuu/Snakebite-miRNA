# vim: set filetype=sh :

rule multiqc_quality_control_mergedData_ars:
    """
    Quality control of fastq files (MULTIQC).
    """
    input:
        unmapped=expand("%s/QC/ars/{samples}_ars_unmapped_fastqc.zip" % (config["project-folder"]), samples=samples),
        mapped=expand("%s/QC/ars/{samples}_ars_mapped_fastqc.zip" % (config["project-folder"]), samples=samples)
    output:
        unmapped=directory("%s/QC/ars/multiqc_unmapped/" % (config["project-folder"])),
        mapped=directory("%s/QC/ars/multiqc_mapped/" % (config["project-folder"]))
    log:
        "%s/logs/MULTIQC/multiqc_ars.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/MULTIQC/multiqc_ars.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        multiqc -f -o {output.unmapped} {input.unmapped} 2> {log};
        multiqc -f -o {output.mapped} {input.mapped} 2> {log};
    """