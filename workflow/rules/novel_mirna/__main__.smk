include: "mpileup.smk"

rule novel_mirna:
    """Predict novel miRNA loci from per-sample genome coverage"""
    input:
        rules.novel_mirna__bedtools__intersect_annotation.output.bed,
        rules.novel_mirna__bedtools__intersect_annotation.output.fasta,
