rule decontaminate__bowtie__phix:
    """Map tRNA-unmapped reads against PhiX to remove sequencing-control contamination (Bowtie)"""
    input:
        reads=rules.decontaminate__bowtie__trna.output.unmapped,
        index=rules.reference__bowtie_index__phix.output,
    output:
        mapped=PHIX_FASTQ / "mapped" / "{sample_id}_PhiX_mapped.fastq",
        unmapped=PHIX_FASTQ / "unmapped" / "{sample_id}_PhiX_unmapped.fastq",
        bam=PHIX_BAM / "{sample_id}_PhiX.bam",
    log:
        "logs/decontaminate/bowtie_phix.{sample_id}.log"
    benchmark:
        "benchmark/decontaminate/bowtie_phix.{sample_id}.tsv"
    params:
        index=features["references"]["phix"],
        m=params["align"]["bowtie"]["m"],
        k=params["align"]["bowtie"]["k"],
    threads: esc("cpus", "decontaminate__bowtie__phix")
    resources:
        runtime=esc("runtime", "decontaminate__bowtie__phix"),
        mem_mb=esc("mem_mb", "decontaminate__bowtie__phix"),
        cpus_per_task=esc("cpus", "decontaminate__bowtie__phix"),
        slurm_partition=esc("partition", "decontaminate__bowtie__phix"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontaminate__bowtie__phix')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontaminate__bowtie__phix"))
    container:
        docker["bowtie"]
    shell:
        """
        exec > {log} 2>&1
        bowtie --best --strata --threads {threads} -k 50 -a -e 99999 --sam \
               --al {output.mapped} --un {output.unmapped} {params.index} {input.reads} \
               | samtools view -bS - > {output.bam}
        touch {output.mapped} {output.unmapped} {output.bam}
        """

rule decontaminate__samtools__phix_flagstat:
    """Report mapping stats for the PhiX-mapped reads (samtools)"""
    input:
        rules.decontaminate__bowtie__phix.output.bam
    output:
        flagstat=PHIX_STATS / "{sample_id}_PhiX.flagstat",
        stats=PHIX_STATS / "{sample_id}_PhiX.stats",
    log:
        "logs/decontaminate/samtools_phix_flagstat.{sample_id}.log"
    benchmark:
        "benchmark/decontaminate/samtools_phix_flagstat.{sample_id}.tsv"
    threads: esc("cpus", "decontaminate__samtools__phix_flagstat")
    resources:
        runtime=esc("runtime", "decontaminate__samtools__phix_flagstat"),
        mem_mb=esc("mem_mb", "decontaminate__samtools__phix_flagstat"),
        cpus_per_task=esc("cpus", "decontaminate__samtools__phix_flagstat"),
        slurm_partition=esc("partition", "decontaminate__samtools__phix_flagstat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontaminate__samtools__phix_flagstat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontaminate__samtools__phix_flagstat"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats}
        """

rule decontaminate__fastqc__phix:
    """Quality control of the PhiX-mapped and PhiX-unmapped reads (FastQC)"""
    input:
        mapped=rules.decontaminate__bowtie__phix.output.mapped,
        unmapped=rules.decontaminate__bowtie__phix.output.unmapped,
    output:
        mapped=QC / "PhiX" / "{sample_id}_PhiX_mapped_fastqc.zip",
        unmapped=QC / "PhiX" / "{sample_id}_PhiX_unmapped_fastqc.zip",
    log:
        QC / "PhiX" / "{sample_id}_fastqc.log"
    benchmark:
        "benchmark/decontaminate/fastqc_phix.{sample_id}.tsv"
    params:
        outfolder=str(QC / "PhiX"),
    threads: esc("cpus", "decontaminate__fastqc__phix")
    resources:
        runtime=esc("runtime", "decontaminate__fastqc__phix"),
        mem_mb=esc("mem_mb", "decontaminate__fastqc__phix"),
        cpus_per_task=esc("cpus", "decontaminate__fastqc__phix"),
        slurm_partition=esc("partition", "decontaminate__fastqc__phix"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontaminate__fastqc__phix')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontaminate__fastqc__phix"))
    container:
        docker["fastqc"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outfolder}
        fastqc -t {threads} -o {params.outfolder} --extract {input.mapped}
        fastqc -t {threads} -o {params.outfolder} --extract {input.unmapped}
        """

rule decontaminate__multiqc__phix:
    """Aggregate the PhiX mapped/unmapped FastQC reports (MultiQC)"""
    input:
        mapped=expand(str(QC / "PhiX" / "{sample_id}_PhiX_mapped_fastqc.zip"), sample_id=SAMPLES),
        unmapped=expand(str(QC / "PhiX" / "{sample_id}_PhiX_unmapped_fastqc.zip"), sample_id=SAMPLES),
    output:
        mapped=directory(QC / "PhiX" / "multiqc_mapped"),
        unmapped=directory(QC / "PhiX" / "multiqc_unmapped"),
    log:
        QC / "PhiX" / "multiqc.log"
    benchmark:
        "benchmark/decontaminate/multiqc_phix.tsv"
    threads: esc("cpus", "decontaminate__multiqc__phix")
    resources:
        runtime=esc("runtime", "decontaminate__multiqc__phix"),
        mem_mb=esc("mem_mb", "decontaminate__multiqc__phix"),
        cpus_per_task=esc("cpus", "decontaminate__multiqc__phix"),
        slurm_partition=esc("partition", "decontaminate__multiqc__phix"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontaminate__multiqc__phix')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontaminate__multiqc__phix"))
    container:
        docker["multiqc"]
    shell:
        """
        exec > {log} 2>&1
        multiqc -f -o {output.mapped} {input.mapped}
        multiqc -f -o {output.unmapped} {input.unmapped}
        """
