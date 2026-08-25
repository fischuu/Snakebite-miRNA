rule report__novel_mirna:
    """Render the novel_mirna module PDF report: novel miRNA locus prediction summary"""
    input:
        rules.novel_mirna.input,
    output:
        pdf=PIPELINE_REPORT / "novel_mirna.pdf",
    log:
        PIPELINE_REPORT / "novel_mirna.log",
    benchmark:
        PIPELINE_REPORT / "novel_mirna_benchmark.tsv",
    container:
        docker["r_report"]
    params:
        script=NOVEL_MIRNA_R,
        pipeline_folder=config["pipeline_folder"],
        features=config["features-file"],
        sample_file=config["sample-file"],
        project_folder=WD,
    threads: esc("cpus", "report__novel_mirna")
    resources:
        runtime=esc("runtime", "report__novel_mirna"),
        mem_mb=esc("mem_mb", "report__novel_mirna"),
        cpus_per_task=esc("cpus", "report__novel_mirna"),
        slurm_partition=esc("partition", "report__novel_mirna"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'report__novel_mirna')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__novel_mirna"))
    shell:
        """
        exec > {log} 2>&1
        R -e "project_folder <- '{params.project_folder}'; \
              pipeline_folder <- '{params.pipeline_folder}'; \
              features_file <- '{params.features}'; \
              sample_file <- '{params.sample_file}'; \
              rmarkdown::render('{params.script}', output_format='pdf_document', \
                                 output_file=file.path('{params.project_folder}', '{output.pdf}'))"
        """
