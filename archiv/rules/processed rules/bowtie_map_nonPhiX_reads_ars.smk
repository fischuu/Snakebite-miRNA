# vim: set filetype=sh :

rule bowtie_map_nonPhiX_reads_ars:
    """
    Map samples that do not map against PhiX against the ARS1.2 bovine genome (bowtie).
    """
    input:
        "%s/FASTQ/PhiX/unmapped/{samples}_PhiX_unmapped.fastq" % (config["project-folder"])
    output:
        mappedReads="%s/FASTQ/ars/mapped/{samples}_ars_mapped.fastq" % (config["project-folder"]),
        unmappedReads="%s/FASTQ/ars/unmapped/{samples}_ars_unmapped.fastq" % (config["project-folder"]),
        file="%s/SAM/ars/{samples}_ars.sam" % (config["project-folder"])
    params:
        index=config["arsIndex"],
        m=config["params"]["bowtie"]["m"],
        k=config["params"]["bowtie"]["k"]
    log:
        "%s/logs/BOWTIE/bowtie_map_ars.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/BOWTIE/bowtie_map_ars.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        bowtie --best --strata -m {params.m} -k {params.k} --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input} > {output.file} 2> {log};
    """
