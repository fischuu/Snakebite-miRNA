# vim: set filetype=sh :

rule bowtie_map_ars_mirna_vs_hg38:
    """
    Map bovine ars miRNA sequences against human hg38 genome(bowtie).
    """
    input:
        "%s/FASTA/ARS/ARS-UCD1.2.95.mirna.fa" % (config["project-folder"])
    output:
        mappedReads="%s/Conservation/mirnaSeq_mapped.fastq" % (config["project-folder"]),
        unmappedReads="%s/Conservation/mirnaSeq_unmapped.fastq" % (config["project-folder"]),
        file="%s/Conservation/mirnaSeq.sam" % (config["project-folder"])
    params:
        index=config["hg38Index"],
        m=config["params"]["bowtie"]["m"],
        k=config["params"]["bowtie"]["k"]
    log:
        "%s/logs/BOWTIE/bowtie_map_ars_mirnaSeq_vs_hg38.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/BOWTIE/bowtie_map_ars_mirnaSeq_vs_hg38.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        bowtie --best --strata -m {params.m} -k {params.k} --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} -f {input} > {output.file} 2> {log};
    """
