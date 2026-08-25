rule report__decontaminate:
    """Render the decontaminate module PDF report: tRNA/PhiX decontamination mapping rates"""
    input:
        rules.decontaminate.input,
    output:
        pdf=PIPELINE_REPORT / "decontaminate.pdf",
    log:
        PIPELINE_REPORT / "decontaminate.log",
    benchmark:
        PIPELINE_REPORT / "decontaminate_benchmark.tsv",
    container:
        docker["r_report"]
    params:
        script=DECONTAMINATE_R,
        pipeline_folder=config["pipeline_folder"],
        features=config["features-file"],
        sample_file=config["sample-file"],
        project_folder=WD,
    threads: esc("cpus", "report__decontaminate")
    resources:
        runtime=esc("runtime", "report__decontaminate"),
        mem_mb=esc("mem_mb", "report__decontaminate"),
        cpus_per_task=esc("cpus", "report__decontaminate"),
        slurm_partition=esc("partition", "report__decontaminate"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'report__decontaminate')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__decontaminate"))
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
