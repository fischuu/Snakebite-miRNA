rule reads__fastqc__run:
    """Quality control of one raw lane's FASTQ (FastQC)"""
    input:
        get_raw_read1
    output:
        READS / "{sample_id}.{library_id}_fastqc.zip",
        html=READS / "{sample_id}.{library_id}_fastqc.html",
    log:
        READS / "{sample_id}.{library_id}_fastqc.log"
    benchmark:
        READS / "benchmark/{sample_id}.{library_id}.tsv"
    params:
        outfolder=str(READS),
        prefix=lambda wc: f"{wc.sample_id}.{wc.library_id}",
    threads: esc("cpus", "reads__fastqc__run")
    resources:
        runtime=esc("runtime", "reads__fastqc__run"),
        mem_mb=esc("mem_mb", "reads__fastqc__run"),
        cpus_per_task=esc("cpus", "reads__fastqc__run"),
        slurm_partition=esc("partition", "reads__fastqc__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reads__fastqc__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reads__fastqc__run"))
    container:
        docker["fastqc"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outfolder}
        fastqc -t {threads} -o {params.outfolder} --extract {input}
        # Normalise FastQC's own naming to the sample_id.library_id prefix
        base=$(basename {input})
        base=${{base%.fastq.gz}}
        base=${{base%.fastq}}
        for ext in zip html; do
            produced="{params.outfolder}/${{base}}_fastqc.${{ext}}"
            wanted="{params.outfolder}/{params.prefix}_fastqc.${{ext}}"
            [ "$produced" != "$wanted" ] && mv "$produced" "$wanted"
        done
        """

rule reads__multiqc__run:
    """Aggregate FastQC reports for every raw lane (MultiQC)"""
    input:
        [
            READS / f"{sample_id}.{library_id}_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
    output:
        directory(READS / "multiqc")
    log:
        READS / "multiqc.log"
    benchmark:
        READS / "benchmark/multiqc.tsv"
    threads: esc("cpus", "reads__multiqc__run")
    resources:
        runtime=esc("runtime", "reads__multiqc__run"),
        mem_mb=esc("mem_mb", "reads__multiqc__run"),
        cpus_per_task=esc("cpus", "reads__multiqc__run"),
        slurm_partition=esc("partition", "reads__multiqc__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reads__multiqc__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reads__multiqc__run"))
    container:
        docker["multiqc"]
    shell:
        """
        exec > {log} 2>&1
        multiqc -f -o {output} {input}
        """
