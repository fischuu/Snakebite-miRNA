include: "__functions__.smk"
include: "bowtie.smk"

rule preparation:
    """Perform the preparations"""
    input:
        rules.preparation__bowtie_trna_index.output,
        
