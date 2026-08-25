include: "mirdb.smk"
include: "bowtie_index.smk"
include: "star_index.smk"

rule reference:
    """Build every Bowtie and STAR index used by the pipeline"""
    input:
        rules.reference__bowtie_index__trna.output,
        rules.reference__bowtie_index__phix.output,
        rules.reference__bowtie_index__mature.output,
        rules.reference__bowtie_index__mature_species.output,
        rules.reference__bowtie_index__hairpin.output,
        rules.reference__bowtie_index__hairpin_species.output,
        rules.reference__bowtie_index__genome.output,
        rules.reference__star_index__mature.output,
        rules.reference__star_index__mature_species.output,
        rules.reference__star_index__hairpin.output,
        rules.reference__star_index__hairpin_species.output,
        rules.reference__star_index__genome.output,
