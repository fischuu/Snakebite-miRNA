# vim: set filetype=sh :

rule fastx_trim_unmapped_ars:
    """
    Trim unmappable ARS sequences to 20 bases (FASTX).
    """
    input:
        "%s/FASTQ/ars/unmapped/{samples}_ars_unmapped.fastq" % (config["project-folder"])
    output:
        "%s/FASTQ/ars_trimmed/{samples}_ars_trimmed.fastq" % (config["project-folder"])
    log:
        "%s/logs/FASTX/fastx.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/FASTX/fastx.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 1
    params:
       length=config["params"]["fastx"]["length"]
    shell:"""
      fastx_trimmer -l {params.length} -i {input} -o {output} 2> {log};
    """