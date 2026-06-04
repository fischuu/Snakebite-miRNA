rule decontamination__samtools_trna_flagstats_run:
    """
    Get mapping stats (samtools).
    """
    input:
        DECON_BOWTIE_TRNA / "bam" / "{samples}.{library}_tRNA.bam"
    output:
        flagstat=DECON_SAMTOOLS_TRNA / "{samples}.{library}_tRNA.flagstat",
        stats=DECON_SAMTOOLS_TRNA / "{samples}.{library}_tRNA.stats",
        sorted=DECON_BOWTIE_TRNA / "bam" / "{samples}.{library}_tRNA.sorted.bam",
        index=DECON_BOWTIE_TRNA / "bam" / "{samples}.{library}_tRNA.sorted.bam.bai",
        idxstats=DECON_SAMTOOLS_TRNA /"{samples}.{library}_tRNA.idxstats"
    log:
        DECON_SAMTOOLS_TRNA / "logs" / "samtools_stats_tRNA.{samples}.{library}.log"
    benchmark:
        DECON_SAMTOOLS_TRNA / "benchmark" / "samtools_stats_tRNA.{samples}.{library}.benchmark.tsv"
    threads: esc("cpus", "decontamination__samtools_trna_flagstats")
    resources:
        runtime=esc("runtime", "decontamination__samtools_trna_flagstats"),
        mem_mb=esc("mem_mb", "decontamination__samtools_trna_flagstats"),
        cpus_per_task=esc("cpus", "decontamination__samtools_trna_flagstats"),
        slurm_partition=esc("partition", "decontamination__samtools_trna_flagstats"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontamination__samtools_trna_flagstats')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontamination__samtools_trna_flagstats"))
    container: docker["samtools"]
    shell:"""
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats};
        
        samtools sort {input} > {output.sorted}
        samtools index {output.sorted}
        samtools idxstats {output.sorted} > {output.idxstats} 
    """


rule samtools_nonTRNA_reads_phix_flagstats:
    """
    Get mapping stats (samtools).
    """
    input:
        DECON_BOWTIE_PHIX / "bam" / "{samples}.{library}_PhiX.bam"
    output:
        flagstat=DECON_SAMTOOLS_PHIX / "{samples}.{library}_PhiX.flagstat",
        stats=DECON_SAMTOOLS_PHIX / "{samples}.{library}_PhiX.stats",
    log:
        DECON_SAMTOOLS_PHIX / "logs" / "samtools_stats_PhiX.{samples}.{library}.log"
    benchmark:
        DECON_SAMTOOLS_PHIX / "benchmark" / "samtools_stats_PhiX.{samples}.{library}.benchmark.tsv"
    threads: esc("cpus", "decontamination__samtools_phix_flagstats")
    resources:
        runtime=esc("runtime", "decontamination__samtools_phix_flagstats"),
        mem_mb=esc("mem_mb", "decontamination__samtools_phix_flagstats"),
        cpus_per_task=esc("cpus", "decontamination__samtools_phix_flagstats"),
        slurm_partition=esc("partition", "decontamination__samtools_phix_flagstats"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontamination__samtools_phix_flagstats')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontamination__samtools_phix_flagstats"))
    container: docker["samtools"]
    shell:"""
        samtools flagstat {input} > {output.flagstat};
        samtools stats {input} > {output.stats};
    """

rule decontamination__samtools_trna_flagstats:
    """Run samtools trna stats"""
    input:
        [
            DECON_SAMTOOLS_TRNA / f"{samples}.{library}_tRNA.flagstat"
            for samples, library in SAMPLES_LIBRARY
        ],
        
rule decontamination__samtools_phix_flagstats:
    """Run samtools phix stats"""
    input:
        [
            DECON_SAMTOOLS_PHIX / f"{samples}.{library}_PhiX.flagstat"
            for samples, library in SAMPLES_LIBRARY
        ],
