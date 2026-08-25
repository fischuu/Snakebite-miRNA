rule align__fastqc__mature:
    """Quality control of the mature Bowtie mapped/unmapped reads (FastQC)"""
    input:
        mapped=rules.align__bowtie__mature.output.mapped,
        unmapped=rules.align__bowtie__mature.output.unmapped,
    output:
        mapped=QC / "Mature" / "{sample_id}_mature_mapped_fastqc.zip",
        unmapped=QC / "Mature" / "{sample_id}_mature_unmapped_fastqc.zip",
    log:
        QC / "Mature" / "{sample_id}_fastqc.log"
    benchmark:
        "benchmark/align/fastqc_mature.{sample_id}.tsv"
    params:
        outfolder=str(QC / "Mature"),
    threads: esc("cpus", "align__fastqc__mature")
    resources:
        runtime=esc("runtime", "align__fastqc__mature"),
        mem_mb=esc("mem_mb", "align__fastqc__mature"),
        cpus_per_task=esc("cpus", "align__fastqc__mature"),
        slurm_partition=esc("partition", "align__fastqc__mature"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__fastqc__mature')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__fastqc__mature"))
    container:
        docker["fastqc"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outfolder}
        fastqc -t {threads} -o {params.outfolder} --extract {input.mapped}
        fastqc -t {threads} -o {params.outfolder} --extract {input.unmapped}
        """

rule align__multiqc__mature:
    """Aggregate the mature mapped/unmapped FastQC reports (MultiQC)"""
    input:
        mapped=expand(str(QC / "Mature" / "{sample_id}_mature_mapped_fastqc.zip"), sample_id=SAMPLES),
        unmapped=expand(str(QC / "Mature" / "{sample_id}_mature_unmapped_fastqc.zip"), sample_id=SAMPLES),
    output:
        mapped=directory(QC / "Mature" / "multiqc_mapped"),
        unmapped=directory(QC / "Mature" / "multiqc_unmapped"),
    log:
        QC / "Mature" / "multiqc.log"
    benchmark:
        "benchmark/align/multiqc_mature.tsv"
    threads: esc("cpus", "align__multiqc__mature")
    resources:
        runtime=esc("runtime", "align__multiqc__mature"),
        mem_mb=esc("mem_mb", "align__multiqc__mature"),
        cpus_per_task=esc("cpus", "align__multiqc__mature"),
        slurm_partition=esc("partition", "align__multiqc__mature"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__multiqc__mature')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__multiqc__mature"))
    container:
        docker["multiqc"]
    shell:
        """
        exec > {log} 2>&1
        multiqc -f -o {output.mapped} {input.mapped}
        multiqc -f -o {output.unmapped} {input.unmapped}
        """

rule align__fastqc__hairpin:
    """Quality control of the hairpin Bowtie mapped/unmapped reads (FastQC)"""
    input:
        mapped=rules.align__bowtie__hairpin.output.mapped,
        unmapped=rules.align__bowtie__hairpin.output.unmapped,
    output:
        mapped=QC / "Hairpin" / "{sample_id}_hairpin_mapped_fastqc.zip",
        unmapped=QC / "Hairpin" / "{sample_id}_hairpin_unmapped_fastqc.zip",
    log:
        QC / "Hairpin" / "{sample_id}_fastqc.log"
    benchmark:
        "benchmark/align/fastqc_hairpin.{sample_id}.tsv"
    params:
        outfolder=str(QC / "Hairpin"),
    threads: esc("cpus", "align__fastqc__hairpin")
    resources:
        runtime=esc("runtime", "align__fastqc__hairpin"),
        mem_mb=esc("mem_mb", "align__fastqc__hairpin"),
        cpus_per_task=esc("cpus", "align__fastqc__hairpin"),
        slurm_partition=esc("partition", "align__fastqc__hairpin"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__fastqc__hairpin')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__fastqc__hairpin"))
    container:
        docker["fastqc"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {params.outfolder}
        fastqc -t {threads} -o {params.outfolder} --extract {input.mapped}
        fastqc -t {threads} -o {params.outfolder} --extract {input.unmapped}
        """

rule align__multiqc__hairpin:
    """Aggregate the hairpin mapped/unmapped FastQC reports (MultiQC)"""
    input:
        mapped=expand(str(QC / "Hairpin" / "{sample_id}_hairpin_mapped_fastqc.zip"), sample_id=SAMPLES),
        unmapped=expand(str(QC / "Hairpin" / "{sample_id}_hairpin_unmapped_fastqc.zip"), sample_id=SAMPLES),
    output:
        mapped=directory(QC / "Hairpin" / "multiqc_mapped"),
        unmapped=directory(QC / "Hairpin" / "multiqc_unmapped"),
    log:
        QC / "Hairpin" / "multiqc.log"
    benchmark:
        "benchmark/align/multiqc_hairpin.tsv"
    threads: esc("cpus", "align__multiqc__hairpin")
    resources:
        runtime=esc("runtime", "align__multiqc__hairpin"),
        mem_mb=esc("mem_mb", "align__multiqc__hairpin"),
        cpus_per_task=esc("cpus", "align__multiqc__hairpin"),
        slurm_partition=esc("partition", "align__multiqc__hairpin"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__multiqc__hairpin')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__multiqc__hairpin"))
    container:
        docker["multiqc"]
    shell:
        """
        exec > {log} 2>&1
        multiqc -f -o {output.mapped} {input.mapped}
        multiqc -f -o {output.unmapped} {input.unmapped}
        """
