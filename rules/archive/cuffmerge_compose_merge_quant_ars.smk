# vim: set filetype=sh :

rule cuffmerge_compose_merge_quant_ars:
    """
    Quantify the mapped reads (cufflinks).
    """
    input:
        expand("%s/GTF/ars/{samples}_ars.gtf" % (config["project-folder"]), samples=samples)
    output:
        txt="%s/filesToMerge_ars.tmp" % (config["project-folder"])
    run:
        with open(output.txt, 'w') as out:
            print(*input, sep="\n", file=out)
