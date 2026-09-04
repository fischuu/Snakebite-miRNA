rule quantify__samtools__mirbase_star:
    """Count STAR primary alignments per miRBase reference (mature/mature_species/hairpin/hairpin_species)"""
    input:
        mature=rules.align__star__mature.output.bam,
        mature_species=rules.align__star__mature_species.output.bam,
        hairpin=rules.align__star__hairpin.output.bam,
        hairpin_species=rules.align__star__hairpin_species.output.bam,
    output:
        mature=QUANT_STAR / "mirbase" / "{sample_id}_star_mature.txt",
        mature_species=QUANT_STAR / "mirbase" / "{sample_id}_star_mature_species.txt",
        hairpin=QUANT_STAR / "mirbase" / "{sample_id}_star_hairpin.txt",
        hairpin_species=QUANT_STAR / "mirbase" / "{sample_id}_star_hairpin_species.txt",
    log:
        QUANT_STAR / "mirbase" / "{sample_id}.log"
    benchmark:
        QUANT_STAR / "mirbase" / "benchmark/{sample_id}.tsv"
    threads: esc("cpus", "quantify__samtools__mirbase_star")
    resources:
        runtime=esc("runtime", "quantify__samtools__mirbase_star"),
        mem_mb=esc("mem_mb", "quantify__samtools__mirbase_star"),
        cpus_per_task=esc("cpus", "quantify__samtools__mirbase_star"),
        slurm_partition=esc("partition", "quantify__samtools__mirbase_star"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'quantify__samtools__mirbase_star')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__samtools__mirbase_star"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools view -F 256 {input.mature} | cut -f3 | sort -T "$(dirname {output.mature})" | uniq -c > {output.mature}
        samtools view -F 256 {input.mature_species} | cut -f3 | sort -T "$(dirname {output.mature_species})" | uniq -c > {output.mature_species}
        samtools view -F 256 {input.hairpin} | cut -f3 | sort -T "$(dirname {output.hairpin})" | uniq -c > {output.hairpin}
        samtools view -F 256 {input.hairpin_species} | cut -f3 | sort -T "$(dirname {output.hairpin_species})" | uniq -c > {output.hairpin_species}
        """

rule quantify__bedtools__novel_mirna:
    """Count reads over the predicted novel-miRNA loci (bedtools multicov)"""
    input:
        bam=rules.align__samtools__remove_secondary.output,
        bed=rules.novel_mirna__bedtools__intersect_annotation.output.bed,
    output:
        QUANT_STAR / "novel_mirna" / "{sample_id}_star_novelMirna_bedtools.txt"
    log:
        QUANT_STAR / "novel_mirna" / "{sample_id}.log"
    benchmark:
        QUANT_STAR / "novel_mirna" / "benchmark/{sample_id}.tsv"
    threads: esc("cpus", "quantify__bedtools__novel_mirna")
    resources:
        runtime=esc("runtime", "quantify__bedtools__novel_mirna"),
        mem_mb=esc("mem_mb", "quantify__bedtools__novel_mirna"),
        cpus_per_task=esc("cpus", "quantify__bedtools__novel_mirna"),
        slurm_partition=esc("partition", "quantify__bedtools__novel_mirna"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'quantify__bedtools__novel_mirna')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__bedtools__novel_mirna"))
    container:
        docker["bedtools"]
    shell:
        """
        exec > {log} 2>&1
        bedtools multicov -bams {input.bam} -bed {input.bed} > {output}
        """
