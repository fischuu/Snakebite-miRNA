rule quantify__seqkit__mature_bowtie:
    """Quantify the Bowtie-aligned reads against mature.fa (seqkit)"""
    input:
        bam=rules.align__bowtie__mature.output.bam,
        bai=rules.align__samtools__mature_flagstat.output.bai,
    output:
        QUANT_BOWTIE / "mature" / "{sample_id}_bowtie_mature_seqkit.txt"
    log:
        QUANT_BOWTIE / "mature" / "{sample_id}.log"
    benchmark:
        QUANT_BOWTIE / "mature" / "benchmark/{sample_id}.tsv"
    threads: esc("cpus", "quantify__seqkit__mature_bowtie")
    resources:
        runtime=esc("runtime", "quantify__seqkit__mature_bowtie"),
        mem_mb=esc("mem_mb", "quantify__seqkit__mature_bowtie"),
        cpus_per_task=esc("cpus", "quantify__seqkit__mature_bowtie"),
        slurm_partition=esc("partition", "quantify__seqkit__mature_bowtie"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'quantify__seqkit__mature_bowtie')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__seqkit__mature_bowtie"))
    container:
        docker["seqkit"]
    shell:
        # seqkit bam -C writes its counts to stderr; capture that instead of treating it as an error
        """
        exec 2> {log}
        seqkit bam -C {input.bam}
        cat {log} > {output}
        """

rule quantify__seqkit__mature_species_bowtie:
    """Quantify the Bowtie-aligned reads against the species-filtered mature.fa (seqkit)"""
    input:
        bam=rules.align__bowtie__mature_species.output.bam,
        bai=rules.align__samtools__mature_species_flagstat.output.bai,
    output:
        QUANT_BOWTIE / "mature_species" / "{sample_id}_bowtie_mature_species_seqkit.txt"
    log:
        QUANT_BOWTIE / "mature_species" / "{sample_id}.log"
    benchmark:
        QUANT_BOWTIE / "mature_species" / "benchmark/{sample_id}.tsv"
    threads: esc("cpus", "quantify__seqkit__mature_species_bowtie")
    resources:
        runtime=esc("runtime", "quantify__seqkit__mature_species_bowtie"),
        mem_mb=esc("mem_mb", "quantify__seqkit__mature_species_bowtie"),
        cpus_per_task=esc("cpus", "quantify__seqkit__mature_species_bowtie"),
        slurm_partition=esc("partition", "quantify__seqkit__mature_species_bowtie"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'quantify__seqkit__mature_species_bowtie')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__seqkit__mature_species_bowtie"))
    container:
        docker["seqkit"]
    shell:
        """
        exec 2> {log}
        seqkit bam -C {input.bam}
        cat {log} > {output}
        """

rule quantify__seqkit__hairpin_bowtie:
    """Quantify the Bowtie-aligned reads against hairpin.fa (seqkit)"""
    input:
        bam=rules.align__bowtie__hairpin.output.bam,
        bai=rules.align__samtools__hairpin_flagstat.output.bai,
    output:
        QUANT_BOWTIE / "hairpin" / "{sample_id}_bowtie_hairpin_seqkit.txt"
    log:
        QUANT_BOWTIE / "hairpin" / "{sample_id}.log"
    benchmark:
        QUANT_BOWTIE / "hairpin" / "benchmark/{sample_id}.tsv"
    threads: esc("cpus", "quantify__seqkit__hairpin_bowtie")
    resources:
        runtime=esc("runtime", "quantify__seqkit__hairpin_bowtie"),
        mem_mb=esc("mem_mb", "quantify__seqkit__hairpin_bowtie"),
        cpus_per_task=esc("cpus", "quantify__seqkit__hairpin_bowtie"),
        slurm_partition=esc("partition", "quantify__seqkit__hairpin_bowtie"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'quantify__seqkit__hairpin_bowtie')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__seqkit__hairpin_bowtie"))
    container:
        docker["seqkit"]
    shell:
        """
        exec 2> {log}
        seqkit bam -C {input.bam}
        cat {log} > {output}
        """
