# vim: set filetype=sh :

rule bowtie_map_unmapped_trimmed_ars_reads_multimapping:
    """
    Map trimmed samples that do not map against ARS against the hg38 human genome (bowtie).
    """
    input:
        "%s/FASTQ/ars_trimmed/unmapped/{samples}_ars_trimmed_unmapped.fastq" % (config["project-folder"])
    output:
        mappedReads="%s/FASTQ/ars_trimmed_multimapped/mapped/{samples}_ars_trimmed_multimapped_mapped.fastq" % (config["project-folder"]),
        unmappedReads="%s/FASTQ/ars_trimmed_multimapped/unmapped/{samples}_ars_trimmed_multimapped_unmapped.fastq" % (config["project-folder"]),
        file="%s/SAM/ars_trimmed_multimapped/{samples}_ars_trimmed_multimapped.sam" % (config["project-folder"])
    params:
        index=config["arsIndex"],
        m=config["params"]["bowtie"]["mMulti"],
        k=config["params"]["bowtie"]["k"]
    log:
        "%s/logs/BOWTIE/bowtie_map_ars_trimmed_multimapped.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/BOWTIE/bowtie_map_ars_trimmed_multimapped.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: config["params"]["bowtie"]["threads"]
    shell:"""
        bowtie --best --strata -m {params.m} -k {params.k} --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input} > {output.file} 2> {log};
    """
