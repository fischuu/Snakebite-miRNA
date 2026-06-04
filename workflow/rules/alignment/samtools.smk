rule alignment__samtools_mature_flagstats_run:
    """
    Get mapping stats (samtools).
    """
    input:
        ALIGN_BOWTIE_MATURE / "bam" / "{samples}.{library}_mature.bam"
    output:
        flagstat=ALIGN_BOWTIE_MATURE / "stats" / "{samples}.{library}_mature.flagstat",
        stats=ALIGN_BOWTIE_MATURE / "stats" / "/{samples}.{library}_mature.stats",
        bai=ALIGN_BOWTIE_MATURE / "bam" / "{samples}.{library}_mature.bam.bai",
    threads: esc("cpus", "alignment__samtools_mature_flagstats")
    resources:
        runtime=esc("runtime", "alignment__samtools_mature_flagstats"),
        mem_mb=esc("mem_mb", "alignment__samtools_mature_flagstats"),
        cpus_per_task=esc("cpus", "alignment__samtools_mature_flagstats"),
        slurm_partition=esc("partition", "alignment__samtools_mature_flagstats"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'alignment__samtools_mature_flagstats')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("alignment__samtools_mature_flagstats"))
    container: docker["samtools"]
    shell:"""
        samtools flagstat {input} > {output.flagstat};
        samtools stats {input} > {output.stats};
        samtools index {input}
    """

rule alignment__samtools_hairpin_flagstats_run:
    """
    Get mapping stats (samtools).
    """
    input:
        ALIGN_BOWTIE_HAIRPIN / "bam" / "{samples}.{library}_hairpin.bam"
    output:
        flagstat=ALIGN_BOWTIE_HAIRPIN / "stats" / "{samples}.{library}_hairpin.flagstat",
        stats=ALIGN_BOWTIE_HAIRPIN / "stats" / "/{samples}.{library}_hairpin.stats",
        bai=ALIGN_BOWTIE_HAIRPIN / "bam" / "{samples}.{library}_hairpin.bam.bai",
    threads: esc("cpus", "alignment__samtools_hairpin_flagstats")
    resources:
        runtime=esc("runtime", "alignment__samtools_hairpin_flagstats"),
        mem_mb=esc("mem_mb", "alignment__samtools_hairpin_flagstats"),
        cpus_per_task=esc("cpus", "alignment__samtools_hairpin_flagstats"),
        slurm_partition=esc("partition", "alignment__samtools_hairpin_flagstats"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'alignment__samtools_hairpin_flagstats')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("alignment__samtools_hairpin_flagstats"))
    container: docker["samtools"]
    shell:"""
        samtools flagstat {input} > {output.flagstat};
        samtools stats {input} > {output.stats};
        samtools index {input}
    """

rule alignment__samtools_mature_species_flagstats_run:
    """
    Get mapping stats (samtools).
    """
    input:
        ALIGN_BOWTIE_MATURE_SPECIES / "bam" / "{samples}.{library}_mature_species.bam"
    output:
        flagstat=ALIGN_BOWTIE_MATURE_SPECIES / "stats" / "{samples}.{library}_mature_species.flagstat",
        stats=ALIGN_BOWTIE_MATURE_SPECIES / "stats" / "/{samples}.{library}_mature_species.stats",
        bai=ALIGN_BOWTIE_MATURE_SPECIES / "bam" / "{samples}.{library}_mature_species.bam.bai",
    threads: esc("cpus", "alignment__samtools_mature_species_flagstats")
    resources:
        runtime=esc("runtime", "alignment__samtools_mature_species_flagstats"),
        mem_mb=esc("mem_mb", "alignment__samtools_mature_species_flagstats"),
        cpus_per_task=esc("cpus", "alignment__samtools_mature_species_flagstats"),
        slurm_partition=esc("partition", "alignment__samtools_mature_species_flagstats"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'alignment__samtools_mature_species_flagstats')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("alignment__samtools_mature_species_flagstats"))
    container: docker["samtools"]
    shell:"""
        samtools flagstat {input} > {output.flagstat};
        samtools stats {input} > {output.stats};
        samtools index {input}
    """

rule alignment__samtools_hairpin_species_flagstats_run:
    """
    Get mapping stats (samtools).
    """
    input:
        ALIGN_BOWTIE_HAIRPIN_SPECIES / "bam" / "{samples}.{library}_hairpin_species.bam"
    output:
        flagstat=ALIGN_BOWTIE_HAIRPIN_SPECIES / "stats" / "{samples}.{library}_hairpin_species.flagstat",
        stats=ALIGN_BOWTIE_HAIRPIN_SPECIES / "stats" / "/{samples}.{library}_hairpin_species.stats",
        bai=ALIGN_BOWTIE_HAIRPIN_SPECIES / "bam" / "{samples}.{library}_hairpin_species.bam.bai",
    threads: esc("cpus", "alignment__samtools_hairpin_species_flagstats")
    resources:
        runtime=esc("runtime", "alignment__samtools_hairpin_species_flagstats"),
        mem_mb=esc("mem_mb", "alignment__samtools_hairpin_species_flagstats"),
        cpus_per_task=esc("cpus", "alignment__samtools_hairpin_species_flagstats"),
        slurm_partition=esc("partition", "alignment__samtools_hairpin_species_flagstats"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'alignment__samtools_hairpin_species_flagstats')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("alignment__samtools_hairpin_species_flagstats"))
    container: docker["samtools"]
    shell:"""
        samtools flagstat {input} > {output.flagstat};
        samtools stats {input} > {output.stats};
        samtools index {input}
    """
    
rule alignment__samtools_reference_flagstats_run:
    """
    Get mapping stats (samtools).
    """
    input:
        ALIGN_BOWTIE_REF / "bam" / "{samples}.{library}_reference.bam"
    output:
        flagstat=ALIGN_BOWTIE_REF / "stats" / "{samples}.{library}_reference.flagstat",
        stats=ALIGN_BOWTIE_REF / "stats" / "/{samples}.{library}_reference.stats",
        bai=ALIGN_BOWTIE_REF / "bam" / "{samples}.{library}_reference.bam.bai",
    threads: esc("cpus", "alignment__samtools_reference_flagstats")
    resources:
        runtime=esc("runtime", "alignment__samtools_reference_flagstats"),
        mem_mb=esc("mem_mb", "alignment__samtools_reference_flagstats"),
        cpus_per_task=esc("cpus", "alignment__samtools_reference_flagstats"),
        slurm_partition=esc("partition", "alignment__samtools_reference_flagstats"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'alignment__samtools_reference_flagstats')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("alignment__samtools_reference_flagstats"))
    container: docker["samtools"]
    shell:"""
        samtools flagstat {input} > {output.flagstat};
        samtools stats {input} > {output.stats};
        samtools index {input}
    """    

rule alignment__samtools_mature_flagstats:
    """Run samtools mature stats"""
    input:
        [
            ALIGN_BOWTIE_MATURE / "stats" / f"{samples}.{library}_mature.flagstat"
            for samples, library in SAMPLES_LIBRARY
        ],

rule alignment__samtools_hairpin_flagstats:
    """Run samtools hairpin stats"""
    input:
        [
            ALIGN_BOWTIE_HAIRPIN / "stats" / f"{samples}.{library}_hairpin.flagstat"
            for samples, library in SAMPLES_LIBRARY
        ],
        
rule alignment__samtools_mature_species_flagstats:
    """Run samtools mature_species stats"""
    input:
        [
            ALIGN_BOWTIE_MATURE_SPECIES / "stats" / f"{samples}.{library}_mature_species.flagstat"
            for samples, library in SAMPLES_LIBRARY
        ],
      
rule alignment__samtools_hairpin_species_flagstats:
    """Run samtools hairpin_species stats"""
    input:
        [
            ALIGN_BOWTIE_HAIRPIN_SPECIES / "stats" / f"{samples}.{library}_hairpin_species.flagstat"
            for samples, library in SAMPLES_LIBRARY
        ],
        
rule alignment__samtools_reference_flagstats:
    """Run samtools reference stats"""
    input:
        [
            ALIGN_BOWTIE_REF / "stats" / f"{samples}.{library}_reference.flagstat"
            for samples, library in SAMPLES_LIBRARY
        ],
