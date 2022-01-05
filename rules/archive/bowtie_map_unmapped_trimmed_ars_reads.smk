# vim: set filetype=sh :

rule bowtie_map_unmapped_trimmed_ars_reads:
    """
    Map trimmed samples that do not map against ARS against the hg38 human genome (bowtie).
    """
    input:
        "%s/FASTQ/ars_trimmed/{samples}_ars_trimmed.fastq" % (config["project-folder"])
    output:
        mappedReads="%s/FASTQ/ars_trimmed/mapped/{samples}_ars_trimmed_mapped.fastq" % (config["project-folder"]),
        unmappedReads="%s/FASTQ/ars_trimmed/unmapped/{samples}_ars_trimmed_unmapped.fastq" % (config["project-folder"]),
        file="%s/SAM/ars_trimmed/{samples}_ars_trimmed.sam" % (config["project-folder"])
    params:
        index=config["arsIndex"],
        m=config["params"]["bowtie"]["m"],
        k=config["params"]["bowtie"]["k"]
    log:
        "%s/logs/BOWTIE/bowtie_map_ars_trimmed.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/BOWTIE/bowtie_map_ars_trimmed.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        bowtie --best --strata -m {params.m} -k {params.k} --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input} > {output.file} 2> {log};
    """
