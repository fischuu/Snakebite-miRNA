include: "__functions__.smk"
include: "cutadapt.smk"

rule preprocess:
    """Perform the preprocessing"""
    input:
        rules.preprocess__cutadapt.input,
        
