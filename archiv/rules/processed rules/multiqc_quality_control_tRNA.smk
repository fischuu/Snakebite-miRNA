# vim: set filetype=sh :

rule multiqc_quality_control_tRNA:
    """
    Quality control of fastq files after tRNA mapping (MULTIQC).
    """
    input:
        unmapped=expand("%s/QC/tRNA/{samples}_tRNA_unmapped_fastqc.zip" % (config["project-folder"]), samples=samples),
        mapped=expand("%s/QC/tRNA/{samples}_tRNA_mapped_fastqc.zip" % (config["project-folder"]), samples=samples)
    output:
        unmapped=directory("%s/QC/tRNA/multiqc_unmapped/" % (config["project-folder"])),
        mapped=directory("%s/QC/tRNA/multiqc_mapped/" % (config["project-folder"]))
    log:
        "%s/logs/MULTIQC/multiqc_tRNA.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/MULTIQC/multiqc_tRNA.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        multiqc -f -o {output.unmapped} {input.unmapped} 2> {log};
        multiqc -f -o {output.mapped} {input.mapped} 2> {log};
    """