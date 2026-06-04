rule preparation__bowtie_trna_index:
    """
    Create tRNA Bowtie Index (BOWTIE).
    """
    input:
        features["references"]["trnaref"]
    output:
        features["references"]["trnaref"] + ".1.ebwt"
    log:
        PREPA_BOWTIE / "log" / "bowtie_index_tRNA.log"
    benchmark:
        PREPA_BOWTIE / "benchmark" / "bowtie_index_tRNA.benchmark.tsv"
    threads: esc("cpus", "preparation__bowtie_trna_index")
    resources:
        runtime=esc("runtime", "preparation__bowtie_trna_index"),
        mem_mb=esc("mem_mb", "preparation__bowtie_trna_index"),
        cpus_per_task=esc("cpus", "preparation__bowtie_trna_index"),
        slurm_partition=esc("partition", "preparation__bowtie_trna_index"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preparation__bowtie_trna_index')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preparation__bowtie_trna_index"))
    container: docker["bowtie"]
    params: 
        output=PREPA_BOWTIE
    shell:"""
     # The input input line is correct, as I want to create the index with the same prefix
      bowtie-build {input} {input} > {log} 2>&1;
    """

    
rule preparation__bowtie_phix_index:
    """
    Create Bowtie Index for phix(BOWTIE).
    """
    input:
        features["references"]["phixref"]
    output:
        features["references"]["phixref"] + ".1.ebwt"
    log:
        PREPA_BOWTIE / "log" / "bowtie_index_PhiX.log"
    benchmark:
        PREPA_BOWTIE / "benchmark" / "bowtie_index_PhiX.benchmark.tsv"
    threads: esc("cpus", "preparation__bowtie_phix_index")
    resources:
        runtime=esc("runtime", "preparation__bowtie_phix_index"),
        mem_mb=esc("mem_mb", "preparation__bowtie_phix_index"),
        cpus_per_task=esc("cpus", "preparation__bowtie_phix_index"),
        slurm_partition=esc("partition", "preparation__bowtie_phix_index"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preparation__bowtie_phix_index')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preparation__bowtie_phix_index"))
    container: docker["bowtie"]
    params: 
        output=PREPA_BOWTIE
    shell:"""
     # The input input line is correct, as I want to create the index with the same prefix
      bowtie-build {input} {input} > {log} 2>&1;
    """

rule preparation__bowtie_mature_index:
    """
    Create Bowtie Index for mature(BOWTIE).
    """
    input:
        REF_MIRBASE / "mature_basesAdjusted.fa",
    output:
        REF_MIRBASE / "mature_basesAdjusted.fa.1.ebwt"
    log:
        REF_MIRBASE / "logs" / "bowtie_index_mature.log"
    benchmark:
        REF_MIRBASE / "benchmark" / "bowtie_index_mature.benchmark.tsv"
    threads: esc("cpus", "preparation__bowtie_mature_index")
    resources:
        runtime=esc("runtime", "preparation__bowtie_mature_index"),
        mem_mb=esc("mem_mb", "preparation__bowtie_mature_index"),
        cpus_per_task=esc("cpus", "preparation__bowtie_mature_index"),
        slurm_partition=esc("partition", "preparation__bowtie_mature_index"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preparation__bowtie_mature_index')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preparation__bowtie_mature_index"))
    container: docker["bowtie"]
    shell:"""
     # The input input line is correct, as I want to create the index with the same prefix
      bowtie-build {input} {input} > {log} 2>&1;
    """

rule preparation__bowtie_hairpin_index:
    """
    Create Bowtie Index for hairpin(BOWTIE).
    """
    input:
        REF_MIRBASE / "hairpin_basesAdjusted.fa",
    output:
        REF_MIRBASE / "hairpin_basesAdjusted.fa.1.ebwt"
    log:
        REF_MIRBASE / "logs" / "bowtie_index_hairpin.log"
    benchmark:
        REF_MIRBASE / "benchmark" / "bowtie_index_hairpin.benchmark.tsv"
    threads: esc("cpus", "preparation__bowtie_hairpin_index")
    resources:
        runtime=esc("runtime", "preparation__bowtie_hairpin_index"),
        mem_mb=esc("mem_mb", "preparation__bowtie_hairpin_index"),
        cpus_per_task=esc("cpus", "preparation__bowtie_hairpin_index"),
        slurm_partition=esc("partition", "preparation__bowtie_hairpin_index"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preparation__bowtie_hairpin_index')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preparation__bowtie_hairpin_index"))
    container: docker["bowtie"]
    shell:"""
     # The input input line is correct, as I want to create the index with the same prefix
      bowtie-build {input} {input} > {log} 2>&1;
    """

rule preparation__bowtie_mature_species_index:
    """
    Create Bowtie Index for mature_species(BOWTIE).
    """
    input:
        REF_MIRBASE / "mature_basesAdjusted_species.fa",
    output:
        REF_MIRBASE / "mature_basesAdjusted_species.fa.1.ebwt"
    log:
        REF_MIRBASE / "logs" / "bowtie_index_mature_species.log"
    benchmark:
        REF_MIRBASE / "benchmark" / "bowtie_index_mature_species.benchmark.tsv"
    threads: esc("cpus", "preparation__bowtie_mature_species_index")
    resources:
        runtime=esc("runtime", "preparation__bowtie_mature_species_index"),
        mem_mb=esc("mem_mb", "preparation__bowtie_mature_species_index"),
        cpus_per_task=esc("cpus", "preparation__bowtie_mature_species_index"),
        slurm_partition=esc("partition", "preparation__bowtie_mature_species_index"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preparation__bowtie_mature_species_index')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preparation__bowtie_mature_species_index"))
    container: docker["bowtie"]
    shell:"""
     # The input input line is correct, as I want to create the index with the same prefix
      bowtie-build {input} {input} > {log} 2>&1;
    """

rule preparation__bowtie_hairpin_species_index:
    """
    Create Bowtie Index for hairpin_species(BOWTIE).
    """
    input:
        REF_MIRBASE / "hairpin_basesAdjusted_species.fa",
    output:
        REF_MIRBASE / "hairpin_basesAdjusted_species.fa.1.ebwt"
    log:
        REF_MIRBASE / "logs" / "bowtie_index_hairpin_species.log"
    benchmark:
        REF_MIRBASE / "benchmark" / "bowtie_index_hairpin_species.benchmark.tsv"
    threads: esc("cpus", "preparation__bowtie_hairpin_species_index")
    resources:
        runtime=esc("runtime", "preparation__bowtie_hairpin_species_index"),
        mem_mb=esc("mem_mb", "preparation__bowtie_hairpin_species_index"),
        cpus_per_task=esc("cpus", "preparation__bowtie_haipin_species_index"),
        slurm_partition=esc("partition", "preparation__bowtie_hairpin_species_index"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preparation__bowtie_hairpin_species_index')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preparation__bowtie_hairpin_species_index"))
    container: docker["bowtie"]
    shell:"""
     # The input input line is correct, as I want to create the index with the same prefix
      bowtie-build {input} {input} > {log} 2>&1;
    """
