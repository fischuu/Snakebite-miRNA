# vim: set filetype=sh :

rule sam_to_sortedBam_ars_trimmed:
    """
    Take the SAM files and prepare the sorted BAMSs (samtools).
    """
    input:
        "%s/SAM/ars_trimmed/{samples}_ars_trimmed.sam" % (config["project-folder"])
    output:
        "%s/BAM/ars_trimmed/{samples}_ars_trimmed.bam" % (config["project-folder"])
    log:
        "%s/logs/SAMTOOLS/samtools_ars_trimmed.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/SAMTOOLS/samtools_ars_trimmed.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        samtools view -b -S {input} | samtools sort > {output} 2> {log};
        samtools index {output};
    """
