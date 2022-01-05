# vim: set filetype=sh :

rule cutadapt_trimming:
    """
    Trimming adapter sequences (CUTADAPT).
    """
    input:
        "%s/FASTQ/RAW/{raw_samples}.fastq.gz" % (config["project-folder"])
    output:
        "%s/FASTQ/TRIMMED/{raw_samples}_trimmed.fastq.gz" % (config["project-folder"])
    log:
        "%s/logs/CUTADAPT/cutadapt.{raw_samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/CUTADAPT/cutadapt.{raw_samples}.benchmark.tsv" % (config["project-folder"])
    threads: config["params"]["cutadapt"]["threads"]
    params:
        adapter=config["adapter"],
        minLength=config["params"]["cutadapt"]["minLength"],
        qualtrim=config["params"]["cutadapt"]["qualtrim"]
    shell:"""
      cutadapt -a file:{params.adapter} \
               --minimum-length {params.minLength} \
               -j {threads} -q {params.qualtrim} --trim-n \
               -o {output} {input} 2> {log};
    """