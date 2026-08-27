include: "seqkit.smk"
include: "featurecounts.smk"
include: "mirbase.smk"

rule quantify:
    """Quantify miRNA and gene-level expression from the Bowtie/STAR alignments"""
    input:
        expand(str(QUANT_BOWTIE / "mature" / "{sample_id}_bowtie_mature_seqkit.txt"), sample_id=SAMPLES),
        expand(str(QUANT_STAR / "genome" / "{sample_id}_star_reference_fc.txt"), sample_id=SAMPLES),
        expand(str(QUANT_STAR / "genome" / "{sample_id}_star_reference_exon_fc.txt"), sample_id=SAMPLES),
        expand(str(QUANT_STAR / "mirbase" / "{sample_id}_star_mature.txt"), sample_id=SAMPLES),
        expand(str(QUANT_STAR / "novel_mirna" / "{sample_id}_star_novelMirna_bedtools.txt"), sample_id=SAMPLES),
