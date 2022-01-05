# vim: set filetype=sh :

rule featureCounts_quantify_ars_STAR:
    """
    Quantify the mapped reads after merging (featureCounts).
    """
    input:
        bam="%s/BAM/STAR/{samples}_star.bam" % (config["project-folder"])
    output:
        file="%s/GTF/STAR/{samples}_star_fc.txt" % (config["project-folder"])
    log:
        "%s/logs/FC/featureCounts_star.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/FC/featureCounts_star.{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        annot=config["arsAnnot"]
    threads: 16
    params: tmpdir=config["scratch-folder"]
    shell:"""
        featureCounts --primary \
                      -p \
                      -T {threads} \
                      -a {params.annot} \
                      -t exon \
                      -g gene_id \
                      -o {output.file} \
                      {input.bam} 2> {log}
    """
