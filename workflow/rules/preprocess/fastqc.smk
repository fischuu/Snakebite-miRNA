rule preprocess__fastqc__trimmed:
    """Quality control of one trimmed lane's FASTQ (FastQC)"""
    input:
        rules.preprocess__cutadapt__run.output.fastq
    output:
        QC / "TRIMMED" / "{sample_id}.{library_id}_trimmed_fastqc.zip"
    log:
        QC / "TRIMMED" / "{sample_id}.{library_id}_fastqc.log"
    benchmark:
        "benchmark/preprocess/fastqc_trimmed.{sample_id}.{library_id}.tsv"
    params:
        outfolder=str(QC / "TRIMMED"),
    threads: esc("cpus", "preprocess__fastqc__trimmed")
    resources:
        runtime=esc("runtime", "preprocess__fastqc__trimmed"),
        mem_mb=esc("mem_mb", "preprocess__fastqc__trimmed"),
        cpus_per_task=esc("cpus", "preprocess__fastqc__trimmed"),
        slurm_partition=esc("partition", "preprocess__fastqc__trimmed"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__fastqc__trimmed')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__fastqc__trimmed"))
    container:
        docker["fastqc"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outfolder}
        fastqc -t {threads} -o {params.outfolder} --extract {input}
        """

rule preprocess__multiqc__trimmed:
    """Aggregate FastQC reports for every trimmed lane (MultiQC)"""
    input:
        [
            QC / "TRIMMED" / f"{sample_id}.{library_id}_trimmed_fastqc.zip"
            for sample_id, library_id in SAMPLE_LIBRARY
        ]
    output:
        directory(QC / "TRIMMED" / "multiqc")
    log:
        QC / "TRIMMED" / "multiqc.log"
    benchmark:
        "benchmark/preprocess/multiqc_trimmed.tsv"
    threads: esc("cpus", "preprocess__multiqc__trimmed")
    resources:
        runtime=esc("runtime", "preprocess__multiqc__trimmed"),
        mem_mb=esc("mem_mb", "preprocess__multiqc__trimmed"),
        cpus_per_task=esc("cpus", "preprocess__multiqc__trimmed"),
        slurm_partition=esc("partition", "preprocess__multiqc__trimmed"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__multiqc__trimmed')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__multiqc__trimmed"))
    container:
        docker["multiqc"]
    shell:
        """
        exec > {log} 2>&1
        multiqc -f -o {output} {input}
        """

rule preprocess__fastqc__concatenated:
    """Quality control of one sample's concatenated FASTQ (FastQC)"""
    input:
        rules.preprocess__concatenate__run.output.fastq
    output:
        QC / "CONCATENATED" / "{sample_id}_R1_fastqc.zip"
    log:
        QC / "CONCATENATED" / "{sample_id}_fastqc.log"
    benchmark:
        "benchmark/preprocess/fastqc_concatenated.{sample_id}.tsv"
    params:
        outfolder=str(QC / "CONCATENATED"),
    threads: esc("cpus", "preprocess__fastqc__concatenated")
    resources:
        runtime=esc("runtime", "preprocess__fastqc__concatenated"),
        mem_mb=esc("mem_mb", "preprocess__fastqc__concatenated"),
        cpus_per_task=esc("cpus", "preprocess__fastqc__concatenated"),
        slurm_partition=esc("partition", "preprocess__fastqc__concatenated"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__fastqc__concatenated')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__fastqc__concatenated"))
    container:
        docker["fastqc"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outfolder}
        fastqc -t {threads} -o {params.outfolder} --extract {input}
        """

rule preprocess__multiqc__concatenated:
    """Aggregate FastQC reports for every sample's concatenated FASTQ (MultiQC)"""
    input:
        expand(str(QC / "CONCATENATED" / "{sample_id}_R1_fastqc.zip"), sample_id=SAMPLES)
    output:
        directory(QC / "CONCATENATED" / "multiqc")
    log:
        QC / "CONCATENATED" / "multiqc.log"
    benchmark:
        "benchmark/preprocess/multiqc_concatenated.tsv"
    threads: esc("cpus", "preprocess__multiqc__concatenated")
    resources:
        runtime=esc("runtime", "preprocess__multiqc__concatenated"),
        mem_mb=esc("mem_mb", "preprocess__multiqc__concatenated"),
        cpus_per_task=esc("cpus", "preprocess__multiqc__concatenated"),
        slurm_partition=esc("partition", "preprocess__multiqc__concatenated"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__multiqc__concatenated')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__multiqc__concatenated"))
    container:
        docker["multiqc"]
    shell:
        """
        exec > {log} 2>&1
        multiqc -f -o {output} {input}
        """
