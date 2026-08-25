include: "seqkit.smk"
include: "featurecounts.smk"
include: "mirbase.smk"

rule quantify:
    """Quantify miRNA and gene-level expression from the Bowtie/STAR alignments"""
    input:
        expand(str(QUANT_BOWTIE / "Mature" / "{sample_id}_bowtie_mature_seqkit.txt"), sample_id=SAMPLES),
        expand(str(QUANT_STAR / "Reference" / "{sample_id}_star_reference_fc.txt"), sample_id=SAMPLES),
        expand(str(QUANT_STAR / "Reference" / "{sample_id}_star_reference_exon_fc.txt"), sample_id=SAMPLES),
        expand(str(QUANT_STAR / "Mirbase" / "{sample_id}_star_mature.txt"), sample_id=SAMPLES),
        expand(str(QUANT_STAR / "Novel_genes" / "{sample_id}_star_novelMirna_bedtools.txt"), sample_id=SAMPLES),
