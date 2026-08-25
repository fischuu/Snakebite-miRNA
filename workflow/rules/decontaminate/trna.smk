rule decontaminate__bowtie__trna:
    """Map concatenated reads against the tRNA reference to remove tRNA contamination (Bowtie)"""
    input:
        reads=rules.preprocess__concatenate__run.output.fastq,
        index=rules.reference__bowtie_index__trna.output,
    output:
        mapped=TRNA_FASTQ / "mapped" / "{sample_id}_tRNA_mapped.fastq",
        unmapped=TRNA_FASTQ / "unmapped" / "{sample_id}_tRNA_unmapped.fastq",
        bam=TRNA_BAM / "{sample_id}_tRNA.bam",
    log:
        "logs/decontaminate/bowtie_trna.{sample_id}.log"
    benchmark:
        "benchmark/decontaminate/bowtie_trna.{sample_id}.tsv"
    params:
        index=features["references"]["tRNA"],
        m=params["align"]["bowtie"]["m"],
        k=params["align"]["bowtie"]["k"],
    threads: esc("cpus", "decontaminate__bowtie__trna")
    resources:
        runtime=esc("runtime", "decontaminate__bowtie__trna"),
        mem_mb=esc("mem_mb", "decontaminate__bowtie__trna"),
        cpus_per_task=esc("cpus", "decontaminate__bowtie__trna"),
        slurm_partition=esc("partition", "decontaminate__bowtie__trna"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontaminate__bowtie__trna')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontaminate__bowtie__trna"))
    container:
        docker["bowtie"]
    shell:
        """
        exec > {log} 2>&1
        bowtie --best --strata --threads {threads} -t -k 50 -a -e 99999 --sam \
               --al {output.mapped} --un {output.unmapped} {params.index} {input.reads} \
               | samtools view -bS - > {output.bam}
        # bowtie/samtools create no files at all on a 0% alignment rate
        touch {output.mapped} {output.unmapped} {output.bam}
        """

rule decontaminate__samtools__trna_flagstat:
    """Sort, index and report mapping stats for the tRNA-mapped reads (samtools)"""
    input:
        rules.decontaminate__bowtie__trna.output.bam
    output:
        flagstat=TRNA_STATS / "{sample_id}_tRNA.flagstat",
        stats=TRNA_STATS / "{sample_id}_tRNA.stats",
        idxstats=TRNA_STATS / "{sample_id}_tRNA.idxstats",
        sorted_bam=TRNA_BAM / "{sample_id}_tRNA.sorted.bam",
    log:
        "logs/decontaminate/samtools_trna_flagstat.{sample_id}.log"
    benchmark:
        "benchmark/decontaminate/samtools_trna_flagstat.{sample_id}.tsv"
    threads: esc("cpus", "decontaminate__samtools__trna_flagstat")
    resources:
        runtime=esc("runtime", "decontaminate__samtools__trna_flagstat"),
        mem_mb=esc("mem_mb", "decontaminate__samtools__trna_flagstat"),
        cpus_per_task=esc("cpus", "decontaminate__samtools__trna_flagstat"),
        slurm_partition=esc("partition", "decontaminate__samtools__trna_flagstat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontaminate__samtools__trna_flagstat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontaminate__samtools__trna_flagstat"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats}
        samtools sort {input} > {output.sorted_bam}
        samtools index {output.sorted_bam}
        samtools idxstats {output.sorted_bam} > {output.idxstats}
        """

rule decontaminate__fastqc__trna:
    """Quality control of the tRNA-mapped and tRNA-unmapped reads (FastQC)"""
    input:
        mapped=rules.decontaminate__bowtie__trna.output.mapped,
        unmapped=rules.decontaminate__bowtie__trna.output.unmapped,
    output:
        mapped=QC / "tRNA" / "{sample_id}_tRNA_mapped_fastqc.zip",
        unmapped=QC / "tRNA" / "{sample_id}_tRNA_unmapped_fastqc.zip",
    log:
        QC / "tRNA" / "{sample_id}_fastqc.log"
    benchmark:
        "benchmark/decontaminate/fastqc_trna.{sample_id}.tsv"
    params:
        outfolder=str(QC / "tRNA"),
    threads: esc("cpus", "decontaminate__fastqc__trna")
    resources:
        runtime=esc("runtime", "decontaminate__fastqc__trna"),
        mem_mb=esc("mem_mb", "decontaminate__fastqc__trna"),
        cpus_per_task=esc("cpus", "decontaminate__fastqc__trna"),
        slurm_partition=esc("partition", "decontaminate__fastqc__trna"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontaminate__fastqc__trna')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontaminate__fastqc__trna"))
    container:
        docker["fastqc"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outfolder}
        fastqc -t {threads} -o {params.outfolder} --extract {input.mapped}
        fastqc -t {threads} -o {params.outfolder} --extract {input.unmapped}
        """

rule decontaminate__multiqc__trna:
    """Aggregate the tRNA mapped/unmapped FastQC reports (MultiQC)"""
    input:
        mapped=expand(str(QC / "tRNA" / "{sample_id}_tRNA_mapped_fastqc.zip"), sample_id=SAMPLES),
        unmapped=expand(str(QC / "tRNA" / "{sample_id}_tRNA_unmapped_fastqc.zip"), sample_id=SAMPLES),
    output:
        mapped=directory(QC / "tRNA" / "multiqc_mapped"),
        unmapped=directory(QC / "tRNA" / "multiqc_unmapped"),
    log:
        QC / "tRNA" / "multiqc.log"
    benchmark:
        "benchmark/decontaminate/multiqc_trna.tsv"
    threads: esc("cpus", "decontaminate__multiqc__trna")
    resources:
        runtime=esc("runtime", "decontaminate__multiqc__trna"),
        mem_mb=esc("mem_mb", "decontaminate__multiqc__trna"),
        cpus_per_task=esc("cpus", "decontaminate__multiqc__trna"),
        slurm_partition=esc("partition", "decontaminate__multiqc__trna"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontaminate__multiqc__trna')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontaminate__multiqc__trna"))
    container:
        docker["multiqc"]
    shell:
        """
        exec > {log} 2>&1
        multiqc -f -o {output.mapped} {input.mapped}
        multiqc -f -o {output.unmapped} {input.unmapped}
        """
