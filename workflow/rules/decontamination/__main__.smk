include: "__functions__.smk"
include: "bowtie.smk"
include: "samtools.smk"

rule decontamination:
    """Perform the decontamination"""
    input:
        rules.decontamination__bowtie_trna.input,
        rules.decontamination__bowtie_phix.input,
        rules.decontamination__samtools_trna_flagstats.input,
        rules.decontamination__samtools_phix_flagstats.input,
