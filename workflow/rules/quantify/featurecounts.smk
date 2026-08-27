rule quantify__featurecounts__genome_bowtie:
    """Quantify the Bowtie-aligned reads against the genome annotation (featureCounts)"""
    input:
        bam=rules.align__bowtie__genome.output.bam
    output:
        QUANT_BOWTIE / "genome" / "{sample_id}_bowtie_reference_fc.txt"
    log:
        QUANT_BOWTIE / "genome" / "{sample_id}.log"
    benchmark:
        QUANT_BOWTIE / "genome" / "benchmark/{sample_id}.tsv"
    params:
        annotation=features["references"]["annotation"],
        g_option=params["quantify"]["featurecounts"]["g_option"],
        t_option=params["quantify"]["featurecounts"]["t_option"],
    threads: esc("cpus", "quantify__featurecounts__genome_bowtie")
    resources:
        runtime=esc("runtime", "quantify__featurecounts__genome_bowtie"),
        mem_mb=esc("mem_mb", "quantify__featurecounts__genome_bowtie"),
        cpus_per_task=esc("cpus", "quantify__featurecounts__genome_bowtie"),
        slurm_partition=esc("partition", "quantify__featurecounts__genome_bowtie"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'quantify__featurecounts__genome_bowtie')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__featurecounts__genome_bowtie"))
    container:
        docker["subread"]
    shell:
        """
        exec > {log} 2>&1
        featureCounts --primary -M -T {threads} \
                      -a {params.annotation} -t {params.t_option} -g {params.g_option} \
                      -o {output} {input.bam}
        """

rule quantify__featurecounts__genome_star:
    """Quantify the STAR-aligned reads against the genome annotation, at gene and exon level (featureCounts)"""
    input:
        bam=rules.align__star__genome.output.bam
    output:
        gene=QUANT_STAR / "genome" / "{sample_id}_star_reference_fc.txt",
        exon=QUANT_STAR / "genome" / "{sample_id}_star_reference_exon_fc.txt",
    log:
        QUANT_STAR / "genome" / "{sample_id}.log"
    benchmark:
        QUANT_STAR / "genome" / "benchmark/{sample_id}.tsv"
    params:
        annotation=features["references"]["annotation"],
        g_option=params["quantify"]["featurecounts"]["g_option"],
        t_option=params["quantify"]["featurecounts"]["t_option"],
        t_option_exon=params["quantify"]["featurecounts"]["t_option_exon"],
    threads: esc("cpus", "quantify__featurecounts__genome_star")
    resources:
        runtime=esc("runtime", "quantify__featurecounts__genome_star"),
        mem_mb=esc("mem_mb", "quantify__featurecounts__genome_star"),
        cpus_per_task=esc("cpus", "quantify__featurecounts__genome_star"),
        slurm_partition=esc("partition", "quantify__featurecounts__genome_star"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'quantify__featurecounts__genome_star')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("quantify__featurecounts__genome_star"))
    container:
        docker["subread"]
    shell:
        """
        exec > {log} 2>&1
        featureCounts --primary -M -T {threads} \
                      -a {params.annotation} -t {params.t_option} -g {params.g_option} \
                      -o {output.gene} {input.bam}
        featureCounts --primary -M -T {threads} \
                      -a {params.annotation} -t {params.t_option_exon} -g {params.g_option} \
                      -o {output.exon} {input.bam}
        """
