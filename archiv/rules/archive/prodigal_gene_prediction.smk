# vim: set filetype=sh :

rule prodigal_gene_prediction:
    """
    Predict genes from assembled reads (PRODIGAL).
    """
    input:
        "%s/MEGAHIT/ars_trimmed_unmapped/ars_trimmed_unmapped_final.contigs.fa" % (config["project-folder"])
    output:
        gtf="%s/PRODIGAL/ars_unmapped/final.contigs.prodigal.gtf" % (config["project-folder"]),
        fa="%s/PRODIGAL/ars_unmapped/final.contigs.prodigal.fa" % (config["project-folder"]),
    log:
        "%s/logs/prodigal_ars_unmapped.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/prodigal_ars_unmapped.benchmark.tsv" % (config["project-folder"])
    threads: 16
    params:
       outfolder=directory("%s/prodigal_ars_unmapped" % (config["project-folder"]))
    shell:"""
       mkdir -p {params.outfolder}

       prodigal -i {input} -o {output.gtf} -a {output.fa} -p meta -f gff
    """
