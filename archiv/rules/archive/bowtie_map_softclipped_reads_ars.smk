# vim: set filetype=sh :

rule bowtie_map_softclipped_reads_ars:
    """
    Map soft-clipped parts from STAR mappings against the ARS1.2 bovine genome (bowtie).
    """
    input:
        "%s/FASTA/STAR/{samples}_softclipped.fasta.15.gz" % (config["project-folder"])
    output:
        mappedReads="%s/FASTQ/ars/mapped/{samples}_ars_mapped_softclipped.fastq" % (config["project-folder"]),
        unmappedReads="%s/FASTQ/ars/unmapped/{samples}_ars_unmapped_softclipped.fastq" % (config["project-folder"]),
        file="%s/SAM/ars/{samples}_ars_softclipped.sam" % (config["project-folder"])
    params:
        index=config["arsIndex"],
        m=config["params"]["bowtie"]["m"],
        k=config["params"]["bowtie"]["k"],
        threads=config["params"]["bowtie"]["threads"]
    log:
        "%s/logs/BOWTIE/bowtie_map_ars_softclipped.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/BOWTIE/bowtie_map_ars_softclipped.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    shell:"""
        bowtie -f --best --strata -m {params.m} -k {params.k} --threads {params.threads}  --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input} > {output.file} 2> {log};
    """
    
rule softclipped_reads_ars_flagstats:
    """
    Get mapping stats (samtools).
    """
    input:
        "%s/SAM/ars/{samples}_ars_softclipped.sam" % (config["project-folder"])
    output:
        "%s/SAM/ars/{samples}_ars_softclipped.flagstat" % (config["project-folder"])
    shell:"""
        samtools flagstat {input} > {output};
    """
