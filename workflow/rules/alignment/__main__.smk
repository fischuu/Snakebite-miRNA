include: "__functions__.smk"
include: "bowtie.smk"

rule alignment:
    """Perform the alignment"""
    input:
        rules.alignment__bowtie_mature.input,
        rules.alignment__bowtie_hairpin.input,
        rules.alignment__bowtie_mature_species.input,
        rules.alignment__bowtie_hairpin_species.input,
