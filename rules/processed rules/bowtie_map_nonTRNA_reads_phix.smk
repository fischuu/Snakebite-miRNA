# vim: set filetype=sh :

rule bowtie_map_nonTRNA_reads_phix:
    """
    Map non-tRNA samples against the PhiX genome (bowtie).
    """
    input:
        "%s/FASTQ/tRNA/unmapped/{samples}_tRNA_unmapped.fastq" % (config["project-folder"])
    output:
        mappedReads="%s/FASTQ/PhiX/mapped/{samples}_PhiX_mapped.fastq" % (config["project-folder"]),
        unmappedReads="%s/FASTQ/PhiX/unmapped/{samples}_PhiX_unmapped.fastq" % (config["project-folder"]),
        file="%s/SAM/PhiX/{samples}_PhiX.sam" % (config["project-folder"])
    params:
        index=config["phixIndex"],
        m=config["params"]["bowtie"]["m"],
        k=config["params"]["bowtie"]["k"]
    log:
        "%s/logs/BOWTIE/bowtie_map_PhiX.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/BOWTIE/bowtie_map_PhiX.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        touch {output.mappedReads}
        touch {output.unmappedReads}
        touch {output.file}

        bowtie --best --strata -m {params.m} -k {params.k} --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input} > {output.file} 2> {log};
    """

rule nonTRNA_reads_phix_flagstats:
    """
    Get mapping stats (samtools).
    """
    input:
        "%s/SAM/PhiX/{samples}_PhiX.sam" % (config["project-folder"])
    output:
        "%s/SAM/PhiX/{samples}_PhiX.flagstat" % (config["project-folder"])
    shell:"""
        samtools flagstat {input} > {output};
    """
