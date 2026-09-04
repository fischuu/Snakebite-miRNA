rule reference__mirdb__prepare:
    """Adjust miRBase mature.fa/hairpin.fa (U->T) and filter them down to one species"""
    input:
        mature=features["references"]["mature"],
        hairpin=features["references"]["hairpin"],
    output:
        mature=MIRBASE / "mature_basesAdjusted.fa",
        mature_species=MIRBASE / "mature_basesAdjusted_species.fa",
        hairpin=MIRBASE / "hairpin_basesAdjusted.fa",
        hairpin_species=MIRBASE / "hairpin_basesAdjusted_species.fa",
        stats=MIRBASE / "mirdb.stats",
    log:
        MIRBASE / "mirdb_prepare.log"
    benchmark:
        MIRBASE / "benchmark/mirdb_prepare.tsv"
    params:
        species=params["species_id"],
        dedupe_script=os.path.join(SCRIPT_FOLDER, "reference", "dedupe_fasta_names.sh"),
    threads: esc("cpus", "reference__mirdb__prepare")
    resources:
        runtime=esc("runtime", "reference__mirdb__prepare"),
        mem_mb=esc("mem_mb", "reference__mirdb__prepare"),
        cpus_per_task=esc("cpus", "reference__mirdb__prepare"),
        slurm_partition=esc("partition", "reference__mirdb__prepare"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'reference__mirdb__prepare')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("reference__mirdb__prepare"))
    container:
        docker["bowtie"]
    shell:
        """
        exec > {log} 2>&1
        sed '/^[^>]/s/U/T/g' {input.mature} | bash {params.dedupe_script} > {output.mature}
        sed '/^[^>]/s/U/T/g' {input.hairpin} | bash {params.dedupe_script} > {output.hairpin}

        grep -A 1 --no-group-separator '{params.species}' {output.mature} | bash {params.dedupe_script} > {output.mature_species}
        grep -A 1 --no-group-separator '{params.species}' {output.hairpin} | bash {params.dedupe_script} > {output.hairpin_species}

        {{
            grep -ce '>' {output.mature}
            grep -ce '>' {output.mature_species}
            grep -ce '>' {output.hairpin}
            grep -ce '>' {output.hairpin_species}
        }} > {output.stats}
        """
