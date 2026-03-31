# vim: set filetype=sh :

rule bash_compose_assemble_ars_unmapped_reads:
    """
    Merge the unmapped reads (BASH).
    """
    input:
        expand("%s/FASTQ/ars/unmapped/{samples}_ars_unmapped.fastq" % (config["project-folder"]), samples=samples)
    output:
        "%s/FASTQ/ars/forMegahit/all_ars_unmapped.fastq" % (config["project-folder"])
    shell:
        "cat {input} > {output}"
