# vim: set filetype=sh :

rule cufflinks_quantify_star:
    """
    Quantify the mapped reads (cufflinks).
    """
    input:
        "%s/BAM/STAR/{samples}_star.bam" % (config["project-folder"])
    output:
        folder=directory("%s/GTF/STAR_cufflinks/{samples}" % (config["project-folder"])),
        file="%s/GTF/STAR_cufflinks/{samples}_star.gtf" % (config["project-folder"])
    params:
        annot=config["arsAnnot"]
    log:
        "%s/logs/CUFFLINKS/cufflinks_star.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/CUFFLINKS/cufflinks_star.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 20
    shell:"""
        cufflinks -p {threads} \
                  -g {params.annot} \
                  -o {output.folder} \
                  -u --overlap-radius 5 --library-type fr-firststrand \
                  {input} 2> {log}

        mv {output.folder}/transcripts.gtf {output.file};
    """

rule cuffmerge_compose_merge_quant_star:
    """
    Quantify the mapped reads (cufflinks).
    """
    input:
        expand("%s/GTF/STAR_cufflinks/{samples}_star.gtf" % (config["project-folder"]), samples=samples)
    output:
        txt="%s/filesToMerge_star.tmp" % (config["project-folder"])
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)

rule cuffmerge_merge_quant_star:
    """
    Merge cufflinks quantifications ARS (cufflinks).
    """
    input:
        "%s/filesToMerge_star.tmp" % (config["project-folder"])
    output:
        folder=directory("%s/References/gtf_merged_out" % (config["project-folder"])),
        file="%s/References/gtf_merged_star.gtf" % (config["project-folder"])
    params:
        annot=config["arsAnnot"]
    log:
        "%s/logs/CUFFMERGE/cuffmerge_ars.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/CUFFMERGE/cuffmerge_ars.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
       cuffmerge -p {threads} -g {params.annot} -o {output.folder} {input} 2> {log};

      # This removes all lines that do not have a gene name
       sed '/gene_name/!d' {output.folder}/merged.gtf > {output.file}
 
      # This is not needed anymore then
      #  mv {output.folder}/merged.gtf {output.file}
    """
