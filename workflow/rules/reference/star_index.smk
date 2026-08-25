rule reference__star_index__mature:
    """Build the STAR index for the species-adjusted mature.fa"""
    input:
        fasta=rules.reference__mirdb__prepare.output.mature
    output:
        STAR_INDEX_MATURE / "chrName.txt"
    log:
        "logs/reference/star_index_mature.log"
    benchmark:
        "benchmark/reference/star_index_mature.tsv"
    params:
        index=str(STAR_INDEX_MATURE),
        tmpdir=str(STAR_TMP / "Mature"),
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
        mkdir -p {params.index}
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
        "logs/reference/star_index_mature_species.log"
    benchmark:
        "benchmark/reference/star_index_mature_species.tsv"
    params:
        index=str(STAR_INDEX_MATURE_SPECIES),
        tmpdir=str(STAR_TMP / "MatureSpecies"),
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
        mkdir -p {params.index}
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
        "logs/reference/star_index_hairpin.log"
    benchmark:
        "benchmark/reference/star_index_hairpin.tsv"
    params:
        index=str(STAR_INDEX_HAIRPIN),
        tmpdir=str(STAR_TMP / "Hairpin"),
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
        mkdir -p {params.index}
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
        "logs/reference/star_index_hairpin_species.log"
    benchmark:
        "benchmark/reference/star_index_hairpin_species.tsv"
    params:
        index=str(STAR_INDEX_HAIRPIN_SPECIES),
        tmpdir=str(STAR_TMP / "HairpinSpecies"),
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
        mkdir -p {params.index}
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
        "logs/reference/star_index_genome.log"
    benchmark:
        "benchmark/reference/star_index_genome.tsv"
    params:
        index=str(STAR_INDEX_GENOME),
        tmpdir=str(STAR_TMP / "Reference"),
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
        mkdir -p {params.index}
        STAR --runThreadN {threads} --outTmpDir {params.tmpdir} --genomeChrBinNbits 15 \
             --limitGenomeGenerateRAM 65000000000 --runMode genomeGenerate \
             --genomeDir {params.index} --genomeFastaFiles {input.fasta}
        """
