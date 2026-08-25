rule reference__bowtie_index__trna:
    """Build the Bowtie index for the tRNA reference"""
    input:
        features["references"]["tRNA"]
    output:
        f"{features['references']['tRNA']}.1.ebwt"
    log:
        "logs/reference/bowtie_index_trna.log"
    benchmark:
        "benchmark/reference/bowtie_index_trna.tsv"
    threads: esc("cpus", "reference__bowtie_index__trna")
    resources:
        runtime=esc("runtime", "reference__bowtie_index__trna"),
        mem_mb=esc("mem_mb", "reference__bowtie_index__trna"),
        cpus_per_task=esc("cpus", "reference__bowtie_index__trna"),
        slurm_partition=esc("partition", "reference__bowtie_index__trna"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__bowtie_index__trna')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__bowtie_index__trna"))
    container:
        docker["bowtie"]
    shell:
        # The input-as-prefix is intentional: the index is built alongside the reference
        "bowtie-build {input} {input} > {log} 2>&1"

rule reference__bowtie_index__phix:
    """Build the Bowtie index for the PhiX reference"""
    input:
        features["references"]["phix"]
    output:
        f"{features['references']['phix']}.1.ebwt"
    log:
        "logs/reference/bowtie_index_phix.log"
    benchmark:
        "benchmark/reference/bowtie_index_phix.tsv"
    threads: esc("cpus", "reference__bowtie_index__phix")
    resources:
        runtime=esc("runtime", "reference__bowtie_index__phix"),
        mem_mb=esc("mem_mb", "reference__bowtie_index__phix"),
        cpus_per_task=esc("cpus", "reference__bowtie_index__phix"),
        slurm_partition=esc("partition", "reference__bowtie_index__phix"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__bowtie_index__phix')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__bowtie_index__phix"))
    container:
        docker["bowtie"]
    shell:
        "bowtie-build {input} {input} > {log} 2>&1"

rule reference__bowtie_index__mature:
    """Build the Bowtie index for the species-adjusted mature.fa"""
    input:
        rules.reference__mirdb__prepare.output.mature
    output:
        f"{REFDIR}/mature_basesAdjusted.fa.1.ebwt"
    log:
        "logs/reference/bowtie_index_mature.log"
    benchmark:
        "benchmark/reference/bowtie_index_mature.tsv"
    threads: esc("cpus", "reference__bowtie_index__mature")
    resources:
        runtime=esc("runtime", "reference__bowtie_index__mature"),
        mem_mb=esc("mem_mb", "reference__bowtie_index__mature"),
        cpus_per_task=esc("cpus", "reference__bowtie_index__mature"),
        slurm_partition=esc("partition", "reference__bowtie_index__mature"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__bowtie_index__mature')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__bowtie_index__mature"))
    container:
        docker["bowtie"]
    shell:
        "bowtie-build {input} {input} > {log} 2>&1"

rule reference__bowtie_index__mature_species:
    """Build the Bowtie index for the species-filtered mature.fa"""
    input:
        rules.reference__mirdb__prepare.output.mature_species
    output:
        f"{REFDIR}/mature_basesAdjusted_species.fa.1.ebwt"
    log:
        "logs/reference/bowtie_index_mature_species.log"
    benchmark:
        "benchmark/reference/bowtie_index_mature_species.tsv"
    threads: esc("cpus", "reference__bowtie_index__mature_species")
    resources:
        runtime=esc("runtime", "reference__bowtie_index__mature_species"),
        mem_mb=esc("mem_mb", "reference__bowtie_index__mature_species"),
        cpus_per_task=esc("cpus", "reference__bowtie_index__mature_species"),
        slurm_partition=esc("partition", "reference__bowtie_index__mature_species"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__bowtie_index__mature_species')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__bowtie_index__mature_species"))
    container:
        docker["bowtie"]
    shell:
        "bowtie-build {input} {input} > {log} 2>&1"

rule reference__bowtie_index__hairpin:
    """Build the Bowtie index for the species-adjusted hairpin.fa"""
    input:
        rules.reference__mirdb__prepare.output.hairpin
    output:
        f"{REFDIR}/hairpin_basesAdjusted.fa.1.ebwt"
    log:
        "logs/reference/bowtie_index_hairpin.log"
    benchmark:
        "benchmark/reference/bowtie_index_hairpin.tsv"
    threads: esc("cpus", "reference__bowtie_index__hairpin")
    resources:
        runtime=esc("runtime", "reference__bowtie_index__hairpin"),
        mem_mb=esc("mem_mb", "reference__bowtie_index__hairpin"),
        cpus_per_task=esc("cpus", "reference__bowtie_index__hairpin"),
        slurm_partition=esc("partition", "reference__bowtie_index__hairpin"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__bowtie_index__hairpin')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__bowtie_index__hairpin"))
    container:
        docker["bowtie"]
    shell:
        "bowtie-build {input} {input} > {log} 2>&1"

rule reference__bowtie_index__hairpin_species:
    """Build the Bowtie index for the species-filtered hairpin.fa"""
    input:
        rules.reference__mirdb__prepare.output.hairpin_species
    output:
        f"{REFDIR}/hairpin_basesAdjusted_species.fa.1.ebwt"
    log:
        "logs/reference/bowtie_index_hairpin_species.log"
    benchmark:
        "benchmark/reference/bowtie_index_hairpin_species.tsv"
    threads: esc("cpus", "reference__bowtie_index__hairpin_species")
    resources:
        runtime=esc("runtime", "reference__bowtie_index__hairpin_species"),
        mem_mb=esc("mem_mb", "reference__bowtie_index__hairpin_species"),
        cpus_per_task=esc("cpus", "reference__bowtie_index__hairpin_species"),
        slurm_partition=esc("partition", "reference__bowtie_index__hairpin_species"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__bowtie_index__hairpin_species')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__bowtie_index__hairpin_species"))
    container:
        docker["bowtie"]
    shell:
        "bowtie-build {input} {input} > {log} 2>&1"

rule reference__bowtie_index__genome:
    """Build the Bowtie index for the host reference genome"""
    input:
        features["references"]["genome"]
    output:
        f"{features['references']['genome']}.1.ebwt"
    log:
        "logs/reference/bowtie_index_genome.log"
    benchmark:
        "benchmark/reference/bowtie_index_genome.tsv"
    threads: esc("cpus", "reference__bowtie_index__genome")
    resources:
        runtime=esc("runtime", "reference__bowtie_index__genome"),
        mem_mb=esc("mem_mb", "reference__bowtie_index__genome"),
        cpus_per_task=esc("cpus", "reference__bowtie_index__genome"),
        slurm_partition=esc("partition", "reference__bowtie_index__genome"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__bowtie_index__genome')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__bowtie_index__genome"))
    container:
        docker["bowtie"]
    shell:
        "bowtie-build --threads {threads} {input} {input} > {log} 2>&1"
