rule novel_mirna__samtools__mpileup:
    """Compute per-sample coverage over the host reference genome (samtools mpileup)"""
    input:
        bam=rules.align__star__genome.output.bam,
        reference=features["references"]["genome"],
    output:
        mpileup=MPILEUP / "{sample_id}.mpileup",
        unique=MPILEUP / "{sample_id}.mpileup.unique",
        bed=MPILEUP / "{sample_id}.bed",
    log:
        "logs/novel_mirna/mpileup.{sample_id}.log"
    benchmark:
        "benchmark/novel_mirna/mpileup.{sample_id}.tsv"
    params:
        unique_script=os.path.join(SCRIPT_FOLDER, "novel_mirna", "uniqueMpileup.sh"),
        bed_script=os.path.join(SCRIPT_FOLDER, "novel_mirna", "mpileupToBed.sh"),
    threads: esc("cpus", "novel_mirna__samtools__mpileup")
    resources:
        runtime=esc("runtime", "novel_mirna__samtools__mpileup"),
        mem_mb=esc("mem_mb", "novel_mirna__samtools__mpileup"),
        cpus_per_task=esc("cpus", "novel_mirna__samtools__mpileup"),
        slurm_partition=esc("partition", "novel_mirna__samtools__mpileup"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'novel_mirna__samtools__mpileup')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("novel_mirna__samtools__mpileup"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools mpileup -d 100 -B -C 50 -f {input.reference} {input.bam} > {output.mpileup}
        bash {params.unique_script} {output.mpileup} > {output.unique}
        bash {params.bed_script} {output.unique} 10 > {output.bed}
        """

rule novel_mirna__bedtools__join_loci:
    """Merge per-sample coverage loci into candidate novel-miRNA regions (bedtools)"""
    input:
        bed=expand(str(MPILEUP / "{sample_id}.bed"), sample_id=SAMPLES)
    output:
        merged=MPILEUP / "merged.bed",
        joined=MPILEUP / "joinedLoci.bed",
    log:
        "logs/novel_mirna/join_loci.log"
    benchmark:
        "benchmark/novel_mirna/join_loci.tsv"
    params:
        join_script=os.path.join(SCRIPT_FOLDER, "novel_mirna", "getJoinedLoci.sh"),
        cover=params["novel_mirna"]["min_cover"],
    threads: esc("cpus", "novel_mirna__bedtools__join_loci")
    resources:
        runtime=esc("runtime", "novel_mirna__bedtools__join_loci"),
        mem_mb=esc("mem_mb", "novel_mirna__bedtools__join_loci"),
        cpus_per_task=esc("cpus", "novel_mirna__bedtools__join_loci"),
        slurm_partition=esc("partition", "novel_mirna__bedtools__join_loci"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'novel_mirna__bedtools__join_loci')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("novel_mirna__bedtools__join_loci"))
    container:
        docker["bedtools"]
    shell:
        """
        exec > {log} 2>&1
        multiIntersectBed -i {input} > {output.merged}
        bash {params.join_script} {output.merged} {params.cover} \
            | bedtools sort | bedtools merge \
            | awk '{{if ($3-$2 > 15) print}}' > {output.joined}
        """

rule novel_mirna__bedtools__intersect_annotation:
    """Remove already-annotated regions from the candidate novel-miRNA loci (bedtools)"""
    input:
        joined=rules.novel_mirna__bedtools__join_loci.output.joined,
        annotation=features["references"]["annotation"],
        fasta=features["references"]["genome"],
    output:
        bed=REFDIR / "novelLoci.bed",
        fasta=REFDIR / "novelLoci.fa",
        unfiltered=temp(REFDIR / "novelLoci_unfiltered.bed"),
    log:
        "logs/novel_mirna/intersect_annotation.log"
    benchmark:
        "benchmark/novel_mirna/intersect_annotation.tsv"
    params:
        filter_script=os.path.join(SCRIPT_FOLDER, "novel_mirna", "filterNovelMirna.sh"),
    threads: esc("cpus", "novel_mirna__bedtools__intersect_annotation")
    resources:
        runtime=esc("runtime", "novel_mirna__bedtools__intersect_annotation"),
        mem_mb=esc("mem_mb", "novel_mirna__bedtools__intersect_annotation"),
        cpus_per_task=esc("cpus", "novel_mirna__bedtools__intersect_annotation"),
        slurm_partition=esc("partition", "novel_mirna__bedtools__intersect_annotation"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'novel_mirna__bedtools__intersect_annotation')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("novel_mirna__bedtools__intersect_annotation"))
    container:
        docker["bedtools"]
    shell:
        """
        exec > {log} 2>&1
        bedtools intersect -v -a {input.joined} -b {input.annotation} > {output.unfiltered}
        bash {params.filter_script} {output.unfiltered} 16 > {output.bed}
        bedtools getfasta -fi {input.fasta} -bed {output.bed} > {output.fasta}
        """
