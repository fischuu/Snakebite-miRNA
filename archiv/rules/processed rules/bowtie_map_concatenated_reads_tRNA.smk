# vim: set filetype=sh :

rule bowtie_map_concatenated_reads_tRNA:
    """
    Map samples against the tRNA database (bowtie).
    """
    input:
        "%s/FASTQ/CONCATENATED/{samples}_R1.fastq.gz" % (config["project-folder"])
    output:
        mappedReads="%s/FASTQ/tRNA/mapped/{samples}_tRNA_mapped.fastq" % (config["project-folder"]),
        unmappedReads="%s/FASTQ/tRNA/unmapped/{samples}_tRNA_unmapped.fastq" % (config["project-folder"]),
        file="%s/SAM/tRNA/{samples}_tRNA.sam" % (config["project-folder"])
    params:
        index=config["tRNAIndex"],
        m=config["params"]["bowtie"]["m"],
        k=config["params"]["bowtie"]["k"]
    log:
        "%s/logs/BOWTIE/bowtie_map_tRNA.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/BOWTIE/bowtie_map_tRNA.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        touch {output.mappedReads}
        touch {output.unmappedReads}
        touch {output.file}
        bowtie --best --strata -m {params.m} -k {params.k} --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input} > {output.file} 2> {log};
    """

rule concatenated_reads_tRNA_flagstats:
    """
    Get mapping stats (samtools).
    """
    input:
        "%s/SAM/tRNA/{samples}_tRNA.sam" % (config["project-folder"])
    output:
        "%s/SAM/tRNA/{samples}_tRNA.flagstat" % (config["project-folder"])
    shell:"""
        samtools flagstat {input} > {output};
    """
