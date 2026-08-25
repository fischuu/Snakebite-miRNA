rule preprocess__cutadapt__run:
    """Trim adapters and low-quality bases from one raw lane (cutadapt)"""
    input:
        get_raw_read1
    output:
        fastq_adapter=TRIMMED / "{sample_id}.{library_id}_trimmed_onlyAdapter.fastq.gz",
        fastq=TRIMMED / "{sample_id}.{library_id}_trimmed.fastq.gz",
    log:
        "logs/preprocess/cutadapt.{sample_id}.{library_id}.log"
    benchmark:
        "benchmark/preprocess/cutadapt.{sample_id}.{library_id}.tsv"
    params:
        fastq_5p_adapter=str(TRIMMED / "{sample_id}.{library_id}_trimmed_only5pAdapter.fastq.gz"),
        adapter5p=params["preprocess"]["cutadapt"]["adapter5p"],
        adapter3p=params["preprocess"]["cutadapt"]["adapter3p"],
        min_length=params["preprocess"]["cutadapt"]["min_length"],
        qualtrim=params["preprocess"]["cutadapt"]["qualtrim"],
        fiveprimetrim=params["preprocess"]["cutadapt"]["fiveprimetrim"],
        threeprimetrim=params["preprocess"]["cutadapt"]["threeprimetrim"],
    threads: esc("cpus", "preprocess__cutadapt__run")
    resources:
        runtime=esc("runtime", "preprocess__cutadapt__run"),
        mem_mb=esc("mem_mb", "preprocess__cutadapt__run"),
        cpus_per_task=esc("cpus", "preprocess__cutadapt__run"),
        slurm_partition=esc("partition", "preprocess__cutadapt__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__cutadapt__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__cutadapt__run"))
    container:
        docker["cutadapt"]
    shell:
        """
        exec > {log} 2>&1
        if [ "{params.adapter5p}" = "" ]; then
            if [ "{params.adapter3p}" = "" ]; then
                echo "NOTE: no 3' or 5' adapter trimming performed"
                cp {input} {output.fastq_adapter}
            else
                echo "NOTE: no 5' adapter trimming performed"
                cutadapt -a {params.adapter3p} -j {threads} -o {output.fastq_adapter} {input}
            fi
        else
            cutadapt -g {params.adapter5p} -j {threads} -o {params.fastq_5p_adapter} {input}
            if [ "{params.adapter3p}" = "" ]; then
                echo "WARNING: no 3' adapter trimming performed"
                cp {params.fastq_5p_adapter} {output.fastq_adapter}
            else
                cutadapt -a {params.adapter3p} -j {threads} -o {output.fastq_adapter} {params.fastq_5p_adapter}
            fi
        fi

        cutadapt --minimum-length {params.min_length} -j {threads} -q {params.qualtrim} --trim-n \
                 --cut {params.fiveprimetrim} --cut -{params.threeprimetrim} \
                 -o {output.fastq} {output.fastq_adapter}
        """
