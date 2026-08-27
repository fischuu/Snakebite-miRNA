rule align__fastqc__mature:
    """Quality control of the mature Bowtie mapped/unmapped reads (FastQC)"""
    input:
        mapped=rules.align__bowtie__mature.output.mapped,
        unmapped=rules.align__bowtie__mature.output.unmapped,
    output:
        mapped=MATURE_BOWTIE / "{sample_id}_mature_mapped_fastqc.zip",
        unmapped=MATURE_BOWTIE / "{sample_id}_mature_unmapped_fastqc.zip",
    log:
        MATURE_BOWTIE / "{sample_id}_fastqc.log"
    benchmark:
        MATURE_BOWTIE / "benchmark/fastqc.{sample_id}.tsv"
    params:
        outfolder=str(MATURE_BOWTIE),
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
        mapped=expand(str(MATURE_BOWTIE / "{sample_id}_mature_mapped_fastqc.zip"), sample_id=SAMPLES),
        unmapped=expand(str(MATURE_BOWTIE / "{sample_id}_mature_unmapped_fastqc.zip"), sample_id=SAMPLES),
    output:
        mapped=directory(MATURE_BOWTIE / "multiqc_mapped"),
        unmapped=directory(MATURE_BOWTIE / "multiqc_unmapped"),
    log:
        MATURE_BOWTIE / "multiqc.log"
    benchmark:
        MATURE_BOWTIE / "benchmark/multiqc.tsv"
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
        mapped=HAIRPIN_BOWTIE / "{sample_id}_hairpin_mapped_fastqc.zip",
        unmapped=HAIRPIN_BOWTIE / "{sample_id}_hairpin_unmapped_fastqc.zip",
    log:
        HAIRPIN_BOWTIE / "{sample_id}_fastqc.log"
    benchmark:
        HAIRPIN_BOWTIE / "benchmark/fastqc.{sample_id}.tsv"
    params:
        outfolder=str(HAIRPIN_BOWTIE),
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
        mapped=expand(str(HAIRPIN_BOWTIE / "{sample_id}_hairpin_mapped_fastqc.zip"), sample_id=SAMPLES),
        unmapped=expand(str(HAIRPIN_BOWTIE / "{sample_id}_hairpin_unmapped_fastqc.zip"), sample_id=SAMPLES),
    output:
        mapped=directory(HAIRPIN_BOWTIE / "multiqc_mapped"),
        unmapped=directory(HAIRPIN_BOWTIE / "multiqc_unmapped"),
    log:
        HAIRPIN_BOWTIE / "multiqc.log"
    benchmark:
        HAIRPIN_BOWTIE / "benchmark/multiqc.tsv"
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
