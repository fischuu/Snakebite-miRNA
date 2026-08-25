rule report__align:
    """Render the align module PDF report: mature/hairpin/genome alignment rates (Bowtie and STAR)"""
    input:
        rules.align.input,
    output:
        pdf=PIPELINE_REPORT / "align.pdf",
    log:
        PIPELINE_REPORT / "align.log",
    benchmark:
        PIPELINE_REPORT / "align_benchmark.tsv",
    container:
        docker["r_report"]
    params:
        script=ALIGN_R,
        pipeline_folder=config["pipeline_folder"],
        features=config["features-file"],
        sample_file=config["sample-file"],
        project_folder=WD,
    threads: esc("cpus", "report__align")
    resources:
        runtime=esc("runtime", "report__align"),
        mem_mb=esc("mem_mb", "report__align"),
        cpus_per_task=esc("cpus", "report__align"),
        slurm_partition=esc("partition", "report__align"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'report__align')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("report__align"))
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
