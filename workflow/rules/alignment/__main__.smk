include: "__functions__.smk"
include: "bowtie.smk"

rule alignment:
    """Perform the alignment"""
    input:
        rules.alignment__bowtie_mature.input,
        rules.alignment__bowtie_hairpin.input,
        rules.alignment__bowtie_mature_species.input,
        rules.alignment__bowtie_hairpin_species.input,
        rules.alignment__bowtie_reference.input,
        rules.alignment__samtools_mature_flagstats.input,
        rules.alignment__samtools_hairpin_flagstats.input,
        rules.alignment__samtools_mature_species_flagstats.input,
        rules.alignment__samtools_hairpin_species_flagstats.input,
