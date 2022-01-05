# vim: set filetype=sh :

rule featureCounts_quantify_merged_ars_trimmed:
    """
    Quantify the mapped reads after merging (featureCounts).
    """
    input:
        bam="%s/BAM/ars_trimmed/{samples}_ars_trimmed.bam" % (config["project-folder"]),
        gtf="%s/References/gtf_merged_ars.gtf" % (config["project-folder"])
    output:
        file="%s/GTF/ars_trimmed_merged_fc/{samples}_ars_trimmed_merged_fc.txt" % (config["project-folder"])
    log:
        "%s/logs/FC/featureCounts_ars_trimmed_merged.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/FC/featureCounts_ars_trimmed_merged.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 16
    params: tmpdir=config["scratch-folder"]
    shell:"""
        TMPDIR=$LOCAL_SCRATCH
        echo "Set the variable TMPDIR to $TMPDIR"
        featureCounts -p \
                      -T {threads} \
                      -a {input.gtf} \
                      -t exon \
                      -g gene_name \
                      -o {output.file} \
                      {input.bam} 2> {log}
    """
