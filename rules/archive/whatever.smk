
####################################################################################################################################
rule bash_get_softclipped_fasta:
    """
    Extract the soft-clipped parts from the alignments.
    """
    input:
        "%s/BAM/STAR/{samples}_star.bam" % (config["project-folder"])
    output:
        "%s/FASTA/STAR/{samples}_softclipped.fasta.gz" % (config["project-folder"])
    log:
        "%s/logs/bash_getSoftclipped.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/bash_getSoftclipped.{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        l=config["params"]["extractSoftClipped"]["length"],
        pipe=config["pipeline-folder"]
    shell:"""
        {params.pipe}/scripts/extractSoftclipped -l {params.l} {input} > {output} 2> {log}
    """    

# The tool for this needs to be on the system path.
# Clone this repository to get it:
# https://github.com/fischuu/SE-MEI

rule bash_filter_softclipped_fasta:
    """
    Filter the soft-clipped parts from the alignments.
    """
    input:
        "%s/FASTA/STAR/{samples}_softclipped.fasta.gz" % (config["project-folder"])
    output:
        "%s/FASTA/STAR/{samples}_softclipped.fasta.15.gz" % (config["project-folder"])
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
        #TMPDIR=$LOCAL_SCRATCH
        echo "Use the variable TMPDIR as "$TMPDIR
        {params.pipe}/scripts/nixshorts_SE_fa {input} 15 2> {log}
    """   
    
rule star_map_softclipped_reads_ars:
    """
    Remap the softclipped reads to the genome (STAR).
    """
    input:
        index=(config["arsSTARIndex"]),
        fastq="%s/FASTA/STAR/{samples}_softclipped.fasta.15.gz" % (config["project-folder"])
    output:
        file="%s/BAM/STAR_softclipped/{samples}_star_softclipped.bam" % (config["project-folder"]),
        dir=directory("%s/BAM/STAR_softclipped/{samples}" % (config["project-folder"]))
    log:
        "%s/logs/star_map_softclipped.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/star_map_softclipped.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 20
    shell:"""
        mkdir -p {output.dir};
        printf \"%s\t%s\t%s\t%s\t%s\n\" {input.index} {input.fastq} {output} {log} {threads}
      	[ ! -d \"{output.dir}\" ] && mkdir {output.dir}
        STAR --genomeDir {input.index} \
            --readFilesIn {input.fastq} \
            --readFilesCommand zcat \
            --outFilterMismatchNoverLmax 0.05 \
            --outFilterMatchNmin 16 \
            --outFilterScoreMinOverLread 0  \
            --outFilterMatchNminOverLread 0 \
            --outFilterMultimapNmax 10 \
            --chimSegmentMin 20 \
            --chimOutType WithinBAM \
            --alignIntronMax 1 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outMultimapperOrder Random \
            --runThreadN {threads} \
            --outFileNamePrefix {wildcards.samples}_ 2> {log};
            
        mv {wildcards.samples}_Aligned.sortedByCoord.out.bam {output.file}
        mv {wildcards.samples}_Log.final.out {wildcards.samples}_Log.progress.out {wildcards.samples}_Log.out {wildcards.samples}_SJ.out.tab {output.dir}
    """
    
