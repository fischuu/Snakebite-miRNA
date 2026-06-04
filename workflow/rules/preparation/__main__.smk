include: "__functions__.smk"
include: "bash.smk"
include: "bowtie.smk"

rule preparation:
    """Perform the preparations"""
    input:
        rules.preparation__bowtie_trna_index.output,
        
