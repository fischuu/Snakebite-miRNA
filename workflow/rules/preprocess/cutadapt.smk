rule preprocess__cutadapt_run:
    """
    Trimming adapter sequences (CUTADAPT).
    """
    input:
        READS / "{sample}.{library}.fq.gz",
    output:
        fastqAdapter = CUTADAPT / "{samples}_trimmed_onlyAdapter.fastq.gz",
        fastq = CUTADAPT / "{samples}_trimmed.fastq.gz",
        wclin = READS / "{samples}.wcl",
        wcloutAdapter = CUTADAPT / "{samples}_trimmed_onlyAdapter.wcl",
        wclout = CUTADAPT / "{samples}_trimmed.wcl",
        wccin = READS / "{samples}.wcc",
        wccoutAdapter = CUTADAPT / "{samples}_trimmed_onlyAdapter.wcc",
        wccout = CUTADAPT / "{samples}_trimmed.wcc"
    log:
        CUTADAPT / "logs/cutadapt.{samples}.log"
    benchmark:
        CUTADAPT / "benchmark/cutadapt.{samples}.benchmark.tsv"
    params:
        fastq5pAdapter = CUTADAPT / "{samples}_trimmed_only5pAdapter.fastq.gz",    
        adapter5p=params["cutadapt"]["adapter5p"],
        adapter3p=params["cutadapt"]["adapter3p"],
        minLength=params["cutadapt"]["minLength"],
        qualtrim=params["cutadapt"]["qualtrim"],
        fiveprimetrim=params["cutadapt"]["fiveprimetrim"],
        threeprimetrim=params["cutadapt"]["threeprimetrim"]
    threads: esc("cpus", "preprocess__cutadapt")
    resources:
        runtime=esc("runtime", "preprocess__cutadapt"),
        mem_mb=esc("mem_mb", "preprocess__cutadapt"),
        cpus_per_task=esc("cpus", "preprocess__cutadapt"),
        slurm_partition=esc("partition", "preprocess__cutadapt"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'preprocess__cutadapt')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("preprocess__cutadapt"))
    singularity: config["singularity"]["cutadapt"]
    shell:"""
      if [ "{params.adapter5p}" = "" ];
      then
          if [ "{params.adapter3p}" = "" ];
          then
             echo "NOTE!!! No 3' or 5' adapter trimming performed!!!" &>> {log};
             cp {input} {output.fastqAdapter}
          else
              echo "NOTE!!! No 5' adapter trimming performed!!!" &>> {log};
              cutadapt -a {params.adapter3p} \
                       -j {threads} \
                       -o {output.fastqAdapter} {input} &>> {log};
          fi
      else
          cutadapt -g {params.adapter5p} \
                   -j {threads} \
                   -o {params.fastq5pAdapter} {input} &>> {log};
                   
          if [ "{params.adapter3p}" = "" ];
          then
              echo "WARNING!!! No 3' adapter trimming performed!!!" &>> {log};
              cp {params.fastq5pAdapter} {output.fastqAdapter}
          else 
              cutadapt -a {params.adapter3p} \
                       -j {threads} \
                       -o {output.fastqAdapter} {params.fastq5pAdapter} &>> {log};
          fi
      fi      

      cutadapt --minimum-length {params.minLength} \
               -j {threads} -q {params.qualtrim} --trim-n \
               --cut {params.fiveprimetrim} --cut -{params.threeprimetrim} \
               -o {output.fastq} {output.fastqAdapter} &>> {log};
               
      zcat {input} | wc -l > {output.wclin}
      zcat {output.fastqAdapter} | wc -l > {output.wcloutAdapter}
      zcat {output.fastq} | wc -l > {output.wclout}
      
      zcat {input} | sed -n '2~4p' | wc -c  > {output.wccin}
      zcat {output.fastqAdapter} | sed -n '2~4p' | wc -c  > {output.wccoutAdapter}
      zcat {output.fastq} | sed -n '2~4p' | wc -c  > {output.wccout}
    """

rule preprocess__cutadapt:
    """Run cutadapt"""
    input:
        [
            READS / f"{samples}_trimmed.fastq.gz"
            for sample in SAMPLES
        ],
