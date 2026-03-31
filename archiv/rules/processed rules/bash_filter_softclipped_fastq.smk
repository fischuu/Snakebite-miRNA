# vim: set filetype=sh :

# The tool for this needs to be on the system path.
# Clone this repository to get it:
# https://github.com/fischuu/SE-MEI

rule bash_filter_softclipped_fastq:
    """
    Filter the soft-clipped parts from the alignments.
    """
    input:
        "%s/FASTA/STAR/{samples}_softclipped.fasta.gz" % (config["project-folder"])
    output:
        file="%s/FASTA/STAR/{samples}_softclipped.fasta.15.gz" % (config["project-folder"])
    log:
        "%s/logs/bash_filterSoftclipped.{samples}.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/bash_filterSoftclipped.{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        folder=directory("%s/FASTQ/STAR/" % (config["project-folder"])),
        tmpfolder=config["scratch-folder"],
        pipe=config["pipeline-folder"]
    shell:"""
    #    TMPDIR="/scratch/project_2001310"
        TMPDIR=$LOCAL_SCRATCH
        echo "Set the variable TMPDIR to $TMPDIR"
        {params.pipe}/scripts/nixshorts_SE_fa {input} 15 2> {log}
    """