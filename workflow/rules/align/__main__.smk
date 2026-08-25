include: "bowtie.smk"
include: "star.smk"
include: "softclipped.smk"
include: "fastqc.smk"

rule align:
    """Align decontaminated reads against mature/hairpin/genome references (Bowtie and STAR)"""
    input:
        expand(str(MATURE_STATS_BOWTIE / "{sample_id}_mature.flagstat"), sample_id=SAMPLES),
        expand(str(HAIRPIN_STATS_BOWTIE / "{sample_id}_hairpin.flagstat"), sample_id=SAMPLES),
        expand(str(MATURE_STATS_STAR / "{sample_id}_mature.flagstat"), sample_id=SAMPLES),
        expand(str(GENOME_STATS_STAR / "{sample_id}_reference.flagstat"), sample_id=SAMPLES),
        expand(str(SOFTCLIPPED / "Reference_softclipped" / "{sample_id}_reference_softclipped.fasta.gz"), sample_id=SAMPLES),
        rules.align__multiqc__mature.output,
        rules.align__multiqc__hairpin.output,
