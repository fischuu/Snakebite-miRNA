include: "cutadapt.smk"
include: "concatenate.smk"
include: "fastqc.smk"

rule preprocess:
    """Trim, concatenate lanes per sample, and QC the results"""
    input:
        expand(str(TRIMMED / "{sample_id}.{library_id}_trimmed.fastq.gz"), zip,
               sample_id=[s for s, l in SAMPLE_LIBRARY], library_id=[l for s, l in SAMPLE_LIBRARY]),
        expand(str(CONCATENATED / "{sample_id}_R1.fastq.gz"), sample_id=SAMPLES),
        rules.preprocess__multiqc__trimmed.output,
        rules.preprocess__multiqc__concatenated.output,
