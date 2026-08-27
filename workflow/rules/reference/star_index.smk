rule reference__star_index__mature:
    """Build the STAR index for the species-adjusted mature.fa"""
    input:
        fasta=rules.reference__mirdb__prepare.output.mature
    output:
        STAR_INDEX_MATURE / "chrName.txt"
    log:
        STAR_INDEX_MATURE / "star_index.log"
    benchmark:
        STAR_INDEX_MATURE / "benchmark.tsv"
    params:
        index=str(STAR_INDEX_MATURE),
        tmpdir=str(STAR_INDEX_TMP / "mature"),
    threads: esc("cpus", "reference__star_index__mature")
    resources:
        runtime=esc("runtime", "reference__star_index__mature"),
        mem_mb=esc("mem_mb", "reference__star_index__mature"),
        cpus_per_task=esc("cpus", "reference__star_index__mature"),
        slurm_partition=esc("partition", "reference__star_index__mature"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__star_index__mature')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__star_index__mature"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        rm -rf {params.tmpdir}
        mkdir -p {params.index} "$(dirname {params.tmpdir})"
        STAR --runThreadN {threads} --outTmpDir {params.tmpdir} --genomeSAindexNbases 6 \
             --limitGenomeGenerateRAM 50000000000 --runMode genomeGenerate \
             --genomeDir {params.index} --genomeFastaFiles {input.fasta}
        """

rule reference__star_index__mature_species:
    """Build the STAR index for the species-filtered mature.fa"""
    input:
        fasta=rules.reference__mirdb__prepare.output.mature_species
    output:
        STAR_INDEX_MATURE_SPECIES / "chrName.txt"
    log:
        STAR_INDEX_MATURE_SPECIES / "star_index.log"
    benchmark:
        STAR_INDEX_MATURE_SPECIES / "benchmark.tsv"
    params:
        index=str(STAR_INDEX_MATURE_SPECIES),
        tmpdir=str(STAR_INDEX_TMP / "mature_species"),
    threads: esc("cpus", "reference__star_index__mature_species")
    resources:
        runtime=esc("runtime", "reference__star_index__mature_species"),
        mem_mb=esc("mem_mb", "reference__star_index__mature_species"),
        cpus_per_task=esc("cpus", "reference__star_index__mature_species"),
        slurm_partition=esc("partition", "reference__star_index__mature_species"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__star_index__mature_species')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__star_index__mature_species"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        rm -rf {params.tmpdir}
        mkdir -p {params.index} "$(dirname {params.tmpdir})"
        STAR --runThreadN {threads} --outTmpDir {params.tmpdir} --genomeSAindexNbases 6 \
             --limitGenomeGenerateRAM 50000000000 --runMode genomeGenerate \
             --genomeDir {params.index} --genomeFastaFiles {input.fasta}
        """

rule reference__star_index__hairpin:
    """Build the STAR index for the species-adjusted hairpin.fa"""
    input:
        fasta=rules.reference__mirdb__prepare.output.hairpin
    output:
        STAR_INDEX_HAIRPIN / "chrName.txt"
    log:
        STAR_INDEX_HAIRPIN / "star_index.log"
    benchmark:
        STAR_INDEX_HAIRPIN / "benchmark.tsv"
    params:
        index=str(STAR_INDEX_HAIRPIN),
        tmpdir=str(STAR_INDEX_TMP / "hairpin"),
    threads: esc("cpus", "reference__star_index__hairpin")
    resources:
        runtime=esc("runtime", "reference__star_index__hairpin"),
        mem_mb=esc("mem_mb", "reference__star_index__hairpin"),
        cpus_per_task=esc("cpus", "reference__star_index__hairpin"),
        slurm_partition=esc("partition", "reference__star_index__hairpin"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__star_index__hairpin')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__star_index__hairpin"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        rm -rf {params.tmpdir}
        mkdir -p {params.index} "$(dirname {params.tmpdir})"
        STAR --runThreadN {threads} --outTmpDir {params.tmpdir} --genomeSAindexNbases 6 \
             --limitGenomeGenerateRAM 50000000000 --runMode genomeGenerate \
             --genomeDir {params.index} --genomeFastaFiles {input.fasta}
        """

rule reference__star_index__hairpin_species:
    """Build the STAR index for the species-filtered hairpin.fa"""
    input:
        fasta=rules.reference__mirdb__prepare.output.hairpin_species
    output:
        STAR_INDEX_HAIRPIN_SPECIES / "chrName.txt"
    log:
        STAR_INDEX_HAIRPIN_SPECIES / "star_index.log"
    benchmark:
        STAR_INDEX_HAIRPIN_SPECIES / "benchmark.tsv"
    params:
        index=str(STAR_INDEX_HAIRPIN_SPECIES),
        tmpdir=str(STAR_INDEX_TMP / "hairpin_species"),
    threads: esc("cpus", "reference__star_index__hairpin_species")
    resources:
        runtime=esc("runtime", "reference__star_index__hairpin_species"),
        mem_mb=esc("mem_mb", "reference__star_index__hairpin_species"),
        cpus_per_task=esc("cpus", "reference__star_index__hairpin_species"),
        slurm_partition=esc("partition", "reference__star_index__hairpin_species"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__star_index__hairpin_species')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__star_index__hairpin_species"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        rm -rf {params.tmpdir}
        mkdir -p {params.index} "$(dirname {params.tmpdir})"
        STAR --runThreadN {threads} --outTmpDir {params.tmpdir} --genomeSAindexNbases 6 \
             --limitGenomeGenerateRAM 50000000000 --runMode genomeGenerate \
             --genomeDir {params.index} --genomeFastaFiles {input.fasta}
        """

rule reference__star_index__genome:
    """Build the STAR index for the host reference genome"""
    input:
        fasta=features["references"]["genome"]
    output:
        STAR_INDEX_GENOME / "chrName.txt"
    log:
        STAR_INDEX_GENOME / "star_index.log"
    benchmark:
        STAR_INDEX_GENOME / "benchmark.tsv"
    params:
        index=str(STAR_INDEX_GENOME),
        tmpdir=str(STAR_INDEX_TMP / "genome"),
    threads: esc("cpus", "reference__star_index__genome")
    resources:
        runtime=esc("runtime", "reference__star_index__genome"),
        mem_mb=esc("mem_mb", "reference__star_index__genome"),
        cpus_per_task=esc("cpus", "reference__star_index__genome"),
        slurm_partition=esc("partition", "reference__star_index__genome"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__star_index__genome')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__star_index__genome"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        rm -rf {params.tmpdir}
        mkdir -p {params.index} "$(dirname {params.tmpdir})"
        STAR --runThreadN {threads} --outTmpDir {params.tmpdir} --genomeChrBinNbits 15 \
             --limitGenomeGenerateRAM 65000000000 --runMode genomeGenerate \
             --genomeDir {params.index} --genomeFastaFiles {input.fasta}
        """
