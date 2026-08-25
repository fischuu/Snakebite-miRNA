include: "reads.smk"
include: "reference.smk"
include: "preprocess.smk"
include: "decontaminate.smk"
include: "align.smk"
include: "quantify.smk"
include: "novel_mirna.smk"

rule report:
    """Generate the report for every module"""
    input:
        rules.report__reads.output,
        rules.report__reference.output,
        rules.report__preprocess.output,
        rules.report__decontaminate.output,
        rules.report__align.output,
        rules.report__quantify.output,
        rules.report__novel_mirna.output,

rule report_reads:
    """Report reads module"""
    input:
        rules.report__reads.output,

rule report_reference:
    """Report reference module"""
    input:
        rules.report__reference.output,

rule report_preprocess:
    """Report preprocess module"""
    input:
        rules.report__preprocess.output,

rule report_decontaminate:
    """Report decontaminate module"""
    input:
        rules.report__decontaminate.output,

rule report_align:
    """Report align module"""
    input:
        rules.report__align.output,

rule report_quantify:
    """Report quantify module"""
    input:
        rules.report__quantify.output,

rule report_novel_mirna:
    """Report novel_mirna module"""
    input:
        rules.report__novel_mirna.output,
