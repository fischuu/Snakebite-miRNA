# vim: set filetype=sh :

rule fastqc_quality_control_tRNA:
    """
    Quality control of fastq files after tRNA mapping (FASTQC).
    """
    input:
        unmapped="%s/FASTQ/tRNA/unmapped/{samples}_tRNA_unmapped.fastq" % (config["project-folder"]),
        mapped="%s/FASTQ/tRNA/mapped/{samples}_tRNA_mapped.fastq" % (config["project-folder"])
    output:
        unmapped="%s/QC/tRNA/{samples}_tRNA_unmapped_fastqc.zip" % (config["project-folder"]),
        mapped="%s/QC/tRNA/{samples}_tRNA_mapped_fastqc.zip" % (config["project-folder"])
    log:
        "%s/logs/FASTQC/fastqc_tRNA.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/FASTQC/fastqc_tRNA.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 12
    params:
        outfolder="%s/QC/tRNA/" % (config["project-folder"])
    shell:"""
        mkdir -p {params.outfolder};
        fastqc -t {threads} -o {params.outfolder} --extract {input.unmapped} 2> {log};
        fastqc -t {threads} -o {params.outfolder} --extract {input.mapped} 2> {log};
    """