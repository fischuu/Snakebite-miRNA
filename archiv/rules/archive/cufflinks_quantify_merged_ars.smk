# vim: set filetype=sh :

rule cufflinks_quantify_merged_ars:
    """
    Quantify the mapped reads after merging (cufflinks).
    """
    input:
        bam="%s/BAM/ars/{samples}_ars.bam" % (config["project-folder"]),
        gtf="%s/References/gtf_merged_ars.gtf" % (config["project-folder"])
    output:
        folder=directory("%s/GTF/ars_merged/{samples}" % (config["project-folder"])),
        file="%s/GTF/ars_merged/{samples}_ars_merged.gtf" % (config["project-folder"])
    log:
        "%s/logs/CUFFLINKS/cufflinks_ars_merged.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/CUFFLINKS/cufflinks_ars_merged.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        cufflinks -p {threads} \
                  -g {input.gtf} \
                  -o {output.folder} \
                  -u --overlap-radius 5 --library-type fr-firststrand \
                  {input} 2> {log}

        mv {output.folder}/transcripts.gtf {output.file};
    """
