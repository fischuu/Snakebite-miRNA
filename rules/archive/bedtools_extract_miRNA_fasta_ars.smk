# vim: set filetype=sh :

rule bedtools_extract_miRNA_fasta_ars:
    """
    Extract nucleotide sequence from all miRNA in ARS1.2 annotation (bedtools).
    """
    input:
        miRNA=config["arsMirnaBed"],
        genome=config["arsRef"]
    output:
        "%s/FASTA/ARS/ARS-UCD1.2.95.mirna.fa" % (config["project-folder"])
    log:
        "%s/logs/bedtools/bedtools_mirna_ars.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/bedtools/bedtools_mirna_ars.benchmark.tsv" % (config["project-folder"])
    shell:"""
       bedtools getfasta -fi {input.genome} -bed {input.miRNA} -name > {output}
    """