# vim: set filetype=sh :

# The tool for this needs to be on the system path.
# Clone this repository to get it:
# https://github.com/fischuu/SE-MEI

rule bash_get_softclipped_fastq:
    """
    Extract the soft-clipped parts from the alignments.
    """
    input:
        "%s/BAM/STAR/{samples}_star.bam" % (config["project-folder"])
    output:
        "%s/FASTA/STAR/{samples}_softclipped.fasta.gz" % (config["project-folder"])
    log:
        "%s/logs/bash_getSoftclipped.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/bash_getSoftclipped.{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        l=config["params"]["extractSoftClipped"]["length"],
        pipe=config["pipeline-folder"]
    shell:"""
        {params.pipe}/scripts/extractSoftclipped -l {params.l} {input} > {output} 2> {log}
    """