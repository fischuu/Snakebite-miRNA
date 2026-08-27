rule preprocess__concatenate__run:
    """Concatenate every trimmed lane belonging to one sample into a single FASTQ"""
    input:
        get_concatenation_inputs
    output:
        fastq=CONCATENATED / "{sample_id}_R1.fastq.gz",
        report=CONCATENATED / "{sample_id}_R1.fastq.gz.report",
    log:
        CONCATENATED / "{sample_id}.log"
    benchmark:
        CONCATENATED / "benchmark/{sample_id}.tsv"
    threads: esc("cpus", "preprocess__concatenate__run")
    resources:
        runtime=esc("runtime", "preprocess__concatenate__run"),
        mem_mb=esc("mem_mb", "preprocess__concatenate__run"),
        cpus_per_task=esc("cpus", "preprocess__concatenate__run"),
        slurm_partition=esc("partition", "preprocess__concatenate__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__concatenate__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__concatenate__run"))
    shell:
        """
        exec > {log} 2>&1
        echo "Concatenating: {input}"
        cat {input} > {output.fastq}
        ls {input} > {output.report}
        """
