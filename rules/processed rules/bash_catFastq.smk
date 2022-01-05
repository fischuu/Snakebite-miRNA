# vim: set filetype=sh :

rule bash_catFastq:
    """
    Concatenate the different lanes into single files (BASH).
    """
    input:
        expand("%s/FASTQ/TRIMMED/{temp}_trimmed.fastq.gz" % (config["project-folder"]), temp=rawsamples)
    output:
        "%s/FASTQ/CONCATENATED/{samples}_R1.fastq.gz" % (config["project-folder"])
    log:
        "%s/logs/BASH/catFastq_{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/BASH/catFastq_{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 1
    params:
       infolder="%s/FASTQ/TRIMMED" % (config["project-folder"]),
       outfolder="%s/FASTQ/CONCATENATED" % (config["project-folder"])
    shell:"""
        mkdir -p {params.outfolder}
        cat {params.infolder}/{wildcards.samples}*_R1_001_trimmed.fastq.gz > {output} 2> {log}
    """