rule decontaminate__bowtie__phix:
    """Map tRNA-unmapped reads against PhiX to remove sequencing-control contamination (Bowtie)"""
    input:
        reads=rules.decontaminate__bowtie__trna.output.unmapped,
        index=rules.reference__bowtie_index__phix.output,
    output:
        mapped=PHIX / "mapped" / "{sample_id}_PhiX_mapped.fastq",
        unmapped=PHIX / "unmapped" / "{sample_id}_PhiX_unmapped.fastq",
        bam=PHIX / "{sample_id}_PhiX.bam",
    log:
        PHIX / "bowtie.{sample_id}.log"
    benchmark:
        PHIX / "benchmark/bowtie.{sample_id}.tsv"
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
        flagstat=PHIX / "{sample_id}_PhiX.flagstat",
        stats=PHIX / "{sample_id}_PhiX.stats",
    log:
        PHIX / "samtools_flagstat.{sample_id}.log"
    benchmark:
        PHIX / "benchmark/samtools_flagstat.{sample_id}.tsv"
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
        mapped=PHIX / "{sample_id}_PhiX_mapped_fastqc.zip",
        unmapped=PHIX / "{sample_id}_PhiX_unmapped_fastqc.zip",
    log:
        PHIX / "{sample_id}_fastqc.log"
    benchmark:
        PHIX / "benchmark/fastqc.{sample_id}.tsv"
    params:
        outfolder=str(PHIX),
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
        mapped=expand(str(PHIX / "{sample_id}_PhiX_mapped_fastqc.zip"), sample_id=SAMPLES),
        unmapped=expand(str(PHIX / "{sample_id}_PhiX_unmapped_fastqc.zip"), sample_id=SAMPLES),
    output:
        mapped=directory(PHIX / "multiqc_mapped"),
        unmapped=directory(PHIX / "multiqc_unmapped"),
    log:
        PHIX / "multiqc.log"
    benchmark:
        PHIX / "benchmark/multiqc.tsv"
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
