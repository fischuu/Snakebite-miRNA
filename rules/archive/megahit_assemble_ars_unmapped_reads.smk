# vim: set filetype=sh :

rule megahit_assemble_ars_unmapped_reads:
    """
    Denovo assemble merged unmapped reads (MEGAHIT).
    """
    input:
        "%s/FASTQ/ars_trimmed/forMegahit/all_ars_trimmed_unmapped.fastq" % (config["project-folder"])
    output:
        temp_dir=directory("%s/temp/megahit_ars_unmapped" % (config["project-folder"])),
        contigs="%s/MEGAHIT/ars_trimmed_unmapped/ars_trimmed_unmapped_final.contigs.fa" % (config["project-folder"])
    log:
        "%s/logs/MEGAHIT/megahit_ars_unmapped.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/MEGAHIT/megahit_ars_unmapped.benchmark.tsv" % (config["project-folder"])
    threads: 16
    shell:"""
        megahit -r {input} -t {threads} -o {output.temp_dir} 2> {log};
        
        mv {output.temp_dir}/final.contigs.fa {output.contigs};
    """
