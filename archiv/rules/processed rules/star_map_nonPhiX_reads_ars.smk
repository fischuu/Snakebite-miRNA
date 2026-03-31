# vim: set filetype=sh :

rule star_map_nonPhiX_reads_ars:
    """
    Map the samples to the genome (STAR).
    """
    input:
        index=(config["arsSTARIndex"]),
        fastq="%s/FASTQ/PhiX/unmapped/{samples}_PhiX_unmapped.fastq" % (config["project-folder"])
    output:
        file="%s/BAM/STAR/{samples}_star.bam" % (config["project-folder"]),
        dir=directory("%s/BAM/STAR/{samples}" % (config["project-folder"]))
    log:
        "%s/logs/star_map.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/star_map.{samples}.benchmark.tsv" % (config["project-folder"])
    threads: 20
    shell:"""
        mkdir -p {output.dir};
        printf \"%s\t%s\t%s\t%s\t%s\n\" {input.index} {input.fastq} {output} {log} {threads}
      	[ ! -d \"{output.dir}\" ] && mkdir {output.dir}
        STAR --genomeDir {input.index} \
            --readFilesIn {input.fastq} \
            --outFilterMismatchNoverLmax 0.01 \
            --outFilterMatchNmin 15 \
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
    
rule featureCounts_quantify_ars_STAR_softclipped_real:
    """
    Quantify the mapped reads after merging (featureCounts).
    """
    input:
        bam="%s/BAM/STAR_softclipped/{samples}_star_softclipped.bam" % (config["project-folder"])
    output:
        file="%s/GTF/STAR_STAR_softclipped/{samples}_softclipped_fc.txt" % (config["project-folder"])
    log:
        "%s/logs/FC/featureCounts_star_star_softclipped.{samples}.log" % (config["project-folder"])
    benchmark:
        "%s/benchmark/FC/featureCounts_star_star_softclipped.{samples}.benchmark.tsv" % (config["project-folder"])
    params:
        annot=config["arsAnnot"]
    threads: 16
    shell:"""
        featureCounts --primary -p \
                      -T {threads} \
                      -a {params.annot} \
                      -t exon \
                      -g gene_id \
                      -o {output.file} \
                      {input.bam} 2> {log}
    """
