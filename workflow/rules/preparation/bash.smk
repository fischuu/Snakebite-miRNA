rule preparation__bash_prepare_mirdb_references:
    """
    Filter and adjust the mirDB references
    """
    input:
        mature=features["references"]["matureref"],
        hairpin=features["references"]["hairpinref"]
    output:
        mature = REF_MIRBASE / "mature_basesAdjusted.fa",
        matureSpecies = REF_MIRBASE / "mature_basesAdjusted_species.fa",
        hairpin = REF_MIRBASE / "hairpin_basesAdjusted.fa",
        hairpinSpecies = REF_MIRBASE / "hairpin_basesAdjusted_species.fa",
        stats = REF_MIRBASE / "mirbase.stats"
    log:
        REF_MIRBASE / "logs" / "prepare_mirdb_references.log"
    benchmark:
        REF_MIRBASE / "benchmark" / "prepare_mirdb_references.benchmark.tsv"
    threads: esc("cpus", "preparation__bash_prepare_mirdb_references")
    resources:
        runtime=esc("runtime", "preparation__bash_prepare_mirdb_references"),
        mem_mb=esc("mem_mb", "preparation__bash_prepare_mirdb_references"),
        cpus_per_task=esc("cpus", "preparation__bash_prepare_mirdb_references"),
        slurm_partition=esc("partition", "preparation__bash_prepare_mirdb_references"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preparation__bash_prepare_mirdb_references')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preparation__bash_prepare_mirdb_references"))    
    params:
        species=params["species-id"]
    shell:""" 
     echo "# Change bases U to T" >> {log}
     sed '/^[^>]/s/U/T/g' {input.mature} > {output.mature} 2>> {log}
     sed '/^[^>]/s/U/T/g' {input.hairpin} > {output.hairpin} 2>> {log}

     echo "# Filter species" >> {log}
     grep -A 1 --no-group-separator '{params.species}' {output.mature}  > {output.matureSpecies} 2>> {log} || true
     grep -A 1 --no-group-separator '{params.species}' {output.hairpin}  > {output.hairpinSpecies} 2>> {log} || true
     
     # Check if the output files are empty and apply fallback
     if [[ ! -s {output.matureSpecies} ]]; then
         echo "Empty matureSpecies file, use the base filtered instead" >> {log}
         cp {output.mature} {output.matureSpecies} 2>> {log}
     fi

     if [[ ! -s {output.hairpinSpecies} ]]; then
         echo "Empty hairpinSpecies file, use the base filtered instead" >> {log}
         cp {output.hairpin} {output.hairpinSpecies} 2>> {log}
     fi


     echo "# Get the line stats" >> {log}
     grep -e '>' {output.mature} | wc -l >> {output.stats} 2>> {log}
     grep -e '>' {output.matureSpecies} | wc -l >> {output.stats} 2>> {log}
     grep -e '>' {output.hairpin} | wc -l >> {output.stats} 2>> {log}
     grep -e '>' {output.hairpinSpecies} | wc -l >> {output.stats} 2>> {log}
    """
