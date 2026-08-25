rule align__extract_softclipped__run:
    """Extract soft-clipped read portions from the STAR alignments (vendored extractSoftclipped binary)"""
    input:
        genome=rules.align__star__genome.output.bam,
        mature=rules.align__star__mature.output.bam,
        mature_species=rules.align__star__mature_species.output.bam,
    output:
        genome=SOFTCLIPPED / "Reference_softclipped" / "{sample_id}_reference_softclipped.fasta.gz",
        mature=SOFTCLIPPED / "Mature_softclipped" / "{sample_id}_mature_softclipped.fasta.gz",
        mature_species=SOFTCLIPPED / "MatureSpecies_softclipped" / "{sample_id}_mature_species_softclipped.fasta.gz",
    log:
        "logs/align/extract_softclipped.{sample_id}.log"
    benchmark:
        "benchmark/align/extract_softclipped.{sample_id}.tsv"
    params:
        length=params["align"]["extract_softclipped"]["length"],
        binary=os.path.join(SCRIPT_FOLDER, "align", "extractSoftclipped"),
    threads: esc("cpus", "align__extract_softclipped__run")
    resources:
        runtime=esc("runtime", "align__extract_softclipped__run"),
        mem_mb=esc("mem_mb", "align__extract_softclipped__run"),
        cpus_per_task=esc("cpus", "align__extract_softclipped__run"),
        slurm_partition=esc("partition", "align__extract_softclipped__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__extract_softclipped__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__extract_softclipped__run"))
    shell:
        """
        exec > {log} 2>&1
        {params.binary} -l {params.length} {input.genome} > {output.genome}
        {params.binary} -l {params.length} {input.mature} > {output.mature}
        {params.binary} -l {params.length} {input.mature_species} > {output.mature_species}
        """
