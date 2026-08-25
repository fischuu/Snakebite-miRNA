include: "trna.smk"
include: "phix.smk"

rule decontaminate:
    """Remove tRNA and PhiX contamination from the concatenated reads"""
    input:
        expand(str(TRNA_FASTQ / "mapped" / "{sample_id}_tRNA_mapped.fastq"), sample_id=SAMPLES),
        expand(str(PHIX_FASTQ / "mapped" / "{sample_id}_PhiX_mapped.fastq"), sample_id=SAMPLES),
        expand(str(TRNA_STATS / "{sample_id}_tRNA.flagstat"), sample_id=SAMPLES),
        expand(str(PHIX_STATS / "{sample_id}_PhiX.flagstat"), sample_id=SAMPLES),
        rules.decontaminate__multiqc__trna.output,
        rules.decontaminate__multiqc__phix.output,
