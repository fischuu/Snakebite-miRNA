# vim: set filetype=sh :

rule bedtools_extract_miRNA_fasta_hg38:
    """
    Extract nucleotide sequence from all miRNA in hg38 annotation (bedtools).
    """
    input:
        miRNA=config["hg38MirnaBed"],
        genome=config["hg38Ref"]
    output:
        "%s/FASTA/hg38/hg38.95.mirna.fa" % (config["project-folder"])
    log:
        "%s/logs/bedtools/bedtools_mirna_hg38.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/bedtools/bedtools_mirna_hg38.benchmark.tsv" % (config["project-folder"])
    shell:"""
       bedtools getfasta -fi {input.genome} -bed {input.miRNA} -name > {output}
    """