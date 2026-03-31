# vim: set filetype=sh :

rule cuffmerge_merge_quant_ars:
    """
    Merge cufflinks quantifications ARS (cufflinks).
    """
    input:
        "%s/filesToMerge_ars.tmp" % (config["project-folder"])
    output:
        folder=directory("%s/References/gtf_merged_out" % (config["project-folder"])),
        file="%s/References/gtf_merged_ars.gtf" % (config["project-folder"])
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
