# vim: set filetype=sh :

rule sam_to_sortedBam_ars:
    """
    Take the SAM files and prepare the sorted BAMSs (samtools).
    """
    input:
        "%s/SAM/ars/{samples}_ars.sam" % (config["project-folder"])
    output:
        "%s/BAM/ars/{samples}_ars.bam" % (config["project-folder"])
    log:
        "%s/logs/SAMTOOLS/samtools_ars.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/SAMTOOLS/samtools_ars.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        samtools view -b -S {input} | samtools sort > {output} 2> {log};
        samtools index {output};
    """
