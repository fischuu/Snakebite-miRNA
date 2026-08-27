rule align__bowtie__mature:
    """Map PhiX-unmapped reads against mature.fa (Bowtie)"""
    input:
        reads=rules.decontaminate__bowtie__phix.output.unmapped,
        index=rules.reference__bowtie_index__mature.output,
    output:
        mapped=MATURE_BOWTIE / "mapped" / "{sample_id}_mature_mapped.fastq",
        unmapped=MATURE_BOWTIE / "unmapped" / "{sample_id}_mature_unmapped.fastq",
        bam=MATURE_BOWTIE / "{sample_id}_mature.bam",
    log:
        MATURE_BOWTIE / "bowtie.{sample_id}.log"
    benchmark:
        MATURE_BOWTIE / "benchmark/bowtie.{sample_id}.tsv"
    params:
        index=str(MIRBASE / "mature_basesAdjusted.fa"),
        m=params["align"]["bowtie"]["m"],
        k=params["align"]["bowtie"]["k"],
    threads: esc("cpus", "align__bowtie__mature")
    resources:
        runtime=esc("runtime", "align__bowtie__mature"),
        mem_mb=esc("mem_mb", "align__bowtie__mature"),
        cpus_per_task=esc("cpus", "align__bowtie__mature"),
        slurm_partition=esc("partition", "align__bowtie__mature"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__bowtie__mature')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__bowtie__mature"))
    container:
        docker["bowtie"]
    shell:
        """
        exec > {log} 2>&1
        bowtie --best --strata --threads {threads} -k 50 -a -e 99999 --sam \
               --al {output.mapped} --un {output.unmapped} {params.index} {input.reads} \
               | samtools view -bS | samtools sort - > {output.bam}
        touch {output.mapped} {output.unmapped} {output.bam}
        """

rule align__samtools__mature_flagstat:
    """Report mapping stats for the mature Bowtie alignment (samtools)"""
    input:
        rules.align__bowtie__mature.output.bam
    output:
        flagstat=MATURE_BOWTIE / "{sample_id}_mature.flagstat",
        stats=MATURE_BOWTIE / "{sample_id}_mature.stats",
        bai=MATURE_BOWTIE / "{sample_id}_mature.bam.bai",
    log:
        MATURE_BOWTIE / "samtools_flagstat.{sample_id}.log"
    benchmark:
        MATURE_BOWTIE / "benchmark/samtools_flagstat.{sample_id}.tsv"
    threads: esc("cpus", "align__samtools__mature_flagstat")
    resources:
        runtime=esc("runtime", "align__samtools__mature_flagstat"),
        mem_mb=esc("mem_mb", "align__samtools__mature_flagstat"),
        cpus_per_task=esc("cpus", "align__samtools__mature_flagstat"),
        slurm_partition=esc("partition", "align__samtools__mature_flagstat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__samtools__mature_flagstat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__samtools__mature_flagstat"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats}
        samtools index {input}
        """

rule align__bowtie__mature_species:
    """Map PhiX-unmapped reads against the species-filtered mature.fa (Bowtie)"""
    input:
        reads=rules.decontaminate__bowtie__phix.output.unmapped,
        index=rules.reference__bowtie_index__mature_species.output,
    output:
        mapped=MATURE_SPECIES_BOWTIE / "mapped" / "{sample_id}_mature_species_mapped.fastq",
        unmapped=MATURE_SPECIES_BOWTIE / "unmapped" / "{sample_id}_mature_species_unmapped.fastq",
        bam=MATURE_SPECIES_BOWTIE / "{sample_id}_mature_species.bam",
    log:
        MATURE_SPECIES_BOWTIE / "bowtie.{sample_id}.log"
    benchmark:
        MATURE_SPECIES_BOWTIE / "benchmark/bowtie.{sample_id}.tsv"
    params:
        index=str(MIRBASE / "mature_basesAdjusted_species.fa"),
        m=params["align"]["bowtie"]["m"],
        k=params["align"]["bowtie"]["k"],
    threads: esc("cpus", "align__bowtie__mature_species")
    resources:
        runtime=esc("runtime", "align__bowtie__mature_species"),
        mem_mb=esc("mem_mb", "align__bowtie__mature_species"),
        cpus_per_task=esc("cpus", "align__bowtie__mature_species"),
        slurm_partition=esc("partition", "align__bowtie__mature_species"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__bowtie__mature_species')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__bowtie__mature_species"))
    container:
        docker["bowtie"]
    shell:
        """
        exec > {log} 2>&1
        bowtie --best --strata --threads {threads} -k 50 -a -e 99999 --sam \
               --al {output.mapped} --un {output.unmapped} {params.index} {input.reads} \
               | samtools view -bS | samtools sort - > {output.bam}
        touch {output.mapped} {output.unmapped} {output.bam}
        """

rule align__samtools__mature_species_flagstat:
    """Report mapping stats for the species-filtered mature Bowtie alignment (samtools)"""
    input:
        rules.align__bowtie__mature_species.output.bam
    output:
        flagstat=MATURE_SPECIES_BOWTIE / "{sample_id}_mature_species.flagstat",
        stats=MATURE_SPECIES_BOWTIE / "{sample_id}_mature_species.stats",
        bai=MATURE_SPECIES_BOWTIE / "{sample_id}_mature_species.bam.bai",
    log:
        MATURE_SPECIES_BOWTIE / "samtools_flagstat.{sample_id}.log"
    benchmark:
        MATURE_SPECIES_BOWTIE / "benchmark/samtools_flagstat.{sample_id}.tsv"
    threads: esc("cpus", "align__samtools__mature_species_flagstat")
    resources:
        runtime=esc("runtime", "align__samtools__mature_species_flagstat"),
        mem_mb=esc("mem_mb", "align__samtools__mature_species_flagstat"),
        cpus_per_task=esc("cpus", "align__samtools__mature_species_flagstat"),
        slurm_partition=esc("partition", "align__samtools__mature_species_flagstat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__samtools__mature_species_flagstat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__samtools__mature_species_flagstat"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats}
        samtools index {input}
        """

rule align__bowtie__hairpin:
    """Map mature-unmapped reads against hairpin.fa (Bowtie)"""
    input:
        reads=rules.align__bowtie__mature.output.unmapped,
        index=rules.reference__bowtie_index__hairpin.output,
    output:
        mapped=HAIRPIN_BOWTIE / "mapped" / "{sample_id}_hairpin_mapped.fastq",
        unmapped=HAIRPIN_BOWTIE / "unmapped" / "{sample_id}_hairpin_unmapped.fastq",
        bam=HAIRPIN_BOWTIE / "{sample_id}_hairpin.bam",
    log:
        HAIRPIN_BOWTIE / "bowtie.{sample_id}.log"
    benchmark:
        HAIRPIN_BOWTIE / "benchmark/bowtie.{sample_id}.tsv"
    params:
        index=str(MIRBASE / "hairpin_basesAdjusted.fa"),
        m=params["align"]["bowtie"]["m"],
        k=params["align"]["bowtie"]["k"],
    threads: esc("cpus", "align__bowtie__hairpin")
    resources:
        runtime=esc("runtime", "align__bowtie__hairpin"),
        mem_mb=esc("mem_mb", "align__bowtie__hairpin"),
        cpus_per_task=esc("cpus", "align__bowtie__hairpin"),
        slurm_partition=esc("partition", "align__bowtie__hairpin"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__bowtie__hairpin')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__bowtie__hairpin"))
    container:
        docker["bowtie"]
    shell:
        """
        exec > {log} 2>&1
        bowtie --best --strata --threads {threads} -k 50 -a -e 99999 --sam \
               --al {output.mapped} --un {output.unmapped} {params.index} {input.reads} \
               | samtools view -bS | samtools sort - > {output.bam}
        touch {output.mapped} {output.unmapped} {output.bam}
        """

rule align__samtools__hairpin_flagstat:
    """Report mapping stats for the hairpin Bowtie alignment (samtools)"""
    input:
        rules.align__bowtie__hairpin.output.bam
    output:
        flagstat=HAIRPIN_BOWTIE / "{sample_id}_hairpin.flagstat",
        stats=HAIRPIN_BOWTIE / "{sample_id}_hairpin.stats",
        bai=HAIRPIN_BOWTIE / "{sample_id}_hairpin.bam.bai",
    log:
        HAIRPIN_BOWTIE / "samtools_flagstat.{sample_id}.log"
    benchmark:
        HAIRPIN_BOWTIE / "benchmark/samtools_flagstat.{sample_id}.tsv"
    threads: esc("cpus", "align__samtools__hairpin_flagstat")
    resources:
        runtime=esc("runtime", "align__samtools__hairpin_flagstat"),
        mem_mb=esc("mem_mb", "align__samtools__hairpin_flagstat"),
        cpus_per_task=esc("cpus", "align__samtools__hairpin_flagstat"),
        slurm_partition=esc("partition", "align__samtools__hairpin_flagstat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__samtools__hairpin_flagstat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__samtools__hairpin_flagstat"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats}
        samtools index {input}
        """

rule align__bowtie__genome:
    """Map hairpin-unmapped reads against the host reference genome (Bowtie).

    Not part of the default 'align' target: kept available for on-demand use,
    mirroring the original pipeline where this branch was reachable but not
    pulled in by 'rule all'.
    """
    input:
        reads=rules.align__bowtie__hairpin.output.unmapped,
        index=rules.reference__bowtie_index__genome.output,
    output:
        mapped=GENOME_BOWTIE / "mapped" / "{sample_id}_reference_mapped.fastq",
        unmapped=GENOME_BOWTIE / "unmapped" / "{sample_id}_reference_unmapped.fastq",
        bam=GENOME_BOWTIE / "{sample_id}_reference.bam",
    log:
        GENOME_BOWTIE / "bowtie.{sample_id}.log"
    benchmark:
        GENOME_BOWTIE / "benchmark/bowtie.{sample_id}.tsv"
    params:
        index=features["references"]["genome"],
        m=params["align"]["bowtie"]["m"],
        k=params["align"]["bowtie"]["k"],
    threads: esc("cpus", "align__bowtie__genome")
    resources:
        runtime=esc("runtime", "align__bowtie__genome"),
        mem_mb=esc("mem_mb", "align__bowtie__genome"),
        cpus_per_task=esc("cpus", "align__bowtie__genome"),
        slurm_partition=esc("partition", "align__bowtie__genome"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__bowtie__genome')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__bowtie__genome"))
    container:
        docker["bowtie"]
    shell:
        """
        exec > {log} 2>&1
        bowtie --best --strata --threads {threads} -k 50 -a -e 99999 --sam \
               --al {output.mapped} --un {output.unmapped} {params.index} {input.reads} \
               | samtools view -bS | samtools sort - > {output.bam}
        touch {output.mapped} {output.unmapped} {output.bam}
        """

rule align__samtools__genome_flagstat:
    """Report mapping stats for the genome Bowtie alignment (samtools). Not part of the default target."""
    input:
        rules.align__bowtie__genome.output.bam
    output:
        flagstat=GENOME_BOWTIE / "{sample_id}_reference.flagstat",
        stats=GENOME_BOWTIE / "{sample_id}_reference.stats",
        bai=GENOME_BOWTIE / "{sample_id}_reference.bam.bai",
    log:
        GENOME_BOWTIE / "samtools_flagstat.{sample_id}.log"
    benchmark:
        GENOME_BOWTIE / "benchmark/samtools_flagstat.{sample_id}.tsv"
    threads: esc("cpus", "align__samtools__genome_flagstat")
    resources:
        runtime=esc("runtime", "align__samtools__genome_flagstat"),
        mem_mb=esc("mem_mb", "align__samtools__genome_flagstat"),
        cpus_per_task=esc("cpus", "align__samtools__genome_flagstat"),
        slurm_partition=esc("partition", "align__samtools__genome_flagstat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__samtools__genome_flagstat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__samtools__genome_flagstat"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats}
        samtools index {input}
        """
