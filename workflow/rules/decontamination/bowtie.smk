rule decontamination__bowtie_trna_run:
    """
    Map samples against the tRNA database (bowtie).
    """
    input:
        reads = CUTADAPT / "{samples}.{library}_trimmed.fastq.gz",
        index = features["references"]["trnaref"] + ".1.ebwt"
    output:
        mappedReads = DECON_BOWTIE_TRNA / "mapped" / "{samples}.{library}_tRNA_mapped.fastq",
        unmappedReads = DECON_BOWTIE_TRNA / "unmapped" / "{samples}.{library}_tRNA_unmapped.fastq",
        file = DECON_BOWTIE_TRNA / "bam" / "{samples}.{library}_tRNA.bam",
        wclmapped = DECON_BOWTIE_TRNA / "mapped" / "{samples}.{library}_tRNA_mapped.wcl",
        wclunmapped = DECON_BOWTIE_TRNA / "unmapped" / "{samples}.{library}_tRNA_unmapped.wcl",
        wccmapped = DECON_BOWTIE_TRNA / "mapped" / "{samples}.{library}_tRNA_mapped.wcc",
        wccunmapped = DECON_BOWTIE_TRNA / "unmapped" / "{samples}.{library}_tRNA_unmapped.wcc"
    params:
        bamFolder=DECON_BOWTIE_TRNA / "bam",
        fastqFolder=DECON_BOWTIE_TRNA,
        index=features["references"]["trnaref"],
        m=params["decontamination"]["bowtie"]["m"],
        k=params["decontamination"]["bowtie"]["k"]
    log:
        DECON_BOWTIE_TRNA / "logs" / "bowtie_map_tRNA.{samples}.{library}.log"
    benchmark:
        DECON_BOWTIE_TRNA / "benchmark" / "bowtie_map_tRNA.{samples}.{library}.benchmark.tsv"
    threads: esc("cpus", "decontamination__bowtie_trna_run")
    resources:
        runtime=esc("runtime", "decontamination__bowtie_trna_run"),
        mem_mb=esc("mem_mb", "decontamination__bowtie_trna_run"),
        cpus_per_task=esc("cpus", "decontamination__bowtie_trna_run"),
        slurm_partition=esc("partition", "decontamination__bowtie_trna_run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontamination__bowtie_trna_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontamination__bowtie_trna_run"))
    container: docker["bowtie"]
    shell:"""
        mkdir -p {params.bamFolder}
        mkdir -p {params.fastqFolder}

        bowtie --best --strata --threads {threads} -t -k 50 -a -e 99999 --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input.reads} | samtools view -bS - > {output.file} 2> {log};
        
        # I need to do this in case that there is a 0.00% alignment rate. In that case no files are created...
        touch {output.mappedReads}
        touch {output.unmappedReads}
        touch {output.file}
        
        wc -l {output.mappedReads} > {output.wclmapped}
        wc -l {output.unmappedReads} > {output.wclunmapped}
        sed -n '2~4p' {output.mappedReads}| wc -c > {output.wccmapped}
        sed -n '2~4p' {output.unmappedReads}| wc -c  > {output.wccunmapped}
    """
    
rule decontamination__bowtie_phix_run:
    """
    Map non-tRNA samples against the PhiX genome (bowtie).
    """
    input:
        reads = DECON_BOWTIE_TRNA / "unmapped" / "{samples}.{library}_tRNA_unmapped.fastq",
        index = features["references"]["phixref"] + ".1.ebwt"
    output:
        mappedReads = DECON_BOWTIE_PHIX / "mapped" / "{samples}.{library}_PhiX_mapped.fastq",
        unmappedReads = DECON_BOWTIE_PHIX / "unmapped" / "{samples}.{library}_PhiX_unmapped.fastq",
        file = DECON_BOWTIE_PHIX / "bam" / "{samples}.{library}_PhiX.bam",
        wclmapped = DECON_BOWTIE_PHIX / "mapped" / "{samples}.{library}_PhiX_mapped.wcl",
        wclunmapped = DECON_BOWTIE_PHIX / "unmapped" / "{samples}.{library}_PhiX_unmapped.wcl",
        wccmapped = DECON_BOWTIE_PHIX / "mapped" / "{samples}.{library}_PhiX_mapped.wcc",
        wccunmapped = DECON_BOWTIE_PHIX / "unmapped" / "{samples}.{library}_PhiX_unmapped.wcc"
    params:
        bamFolder=DECON_BOWTIE_PHIX / "bam",
        fastqFolder=DECON_BOWTIE_PHIX,
        index=features["references"]["phixref"],
        m=params["decontamination"]["bowtie"]["m"],
        k=params["decontamination"]["bowtie"]["k"]
    log:
        DECON_BOWTIE_PHIX / "logs" / "bowtie_map_PhiX.{samples}.{library}.log"
    benchmark:
        DECON_BOWTIE_PHIX / "benchmark" / "bowtie_map_PhiX.{samples}.{library}.benchmark.tsv"
    threads: esc("cpus", "decontamination__bowtie_phix_run")
    resources:
        runtime=esc("runtime", "decontamination__bowtie_phix_run"),
        mem_mb=esc("mem_mb", "decontamination__bowtie_phix_run"),
        cpus_per_task=esc("cpus", "decontamination__bowtie_phix_run"),
        slurm_partition=esc("partition", "decontamination__bowtie_phix_run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'decontamination__bowtie_phix_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("decontamination__bowtie_phix_run"))
    container: docker["bowtie"]

    shell:"""
        mkdir -p {params.bamFolder}
        mkdir -p {params.fastqFolder}

        bowtie --best --strata --threads {threads} -k 50 -a -e 99999 --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input.reads} | samtools view -bS - > {output.file} 2> {log};
        
        # I need to do this in case that there is a 0.00% alignment rate. In that case no files are created...
        touch {output.mappedReads}
        touch {output.unmappedReads}
        touch {output.file}
        
        wc -l {output.mappedReads} > {output.wclmapped}
        wc -l {output.unmappedReads} > {output.wclunmapped}
        
        sed -n '2~4p' {output.mappedReads}| wc -c > {output.wccmapped}
        sed -n '2~4p' {output.unmappedReads}| wc -c  > {output.wccunmapped}
    """
    

rule decontamination__bowtie_trna:
    """Run bowtie trna"""
    input:
        [
            DECON_BOWTIE_TRNA / "bam" / f"{samples}.{library}_tRNA.bam"
            for samples, library in SAMPLES_LIBRARY
        ],
        
rule decontamination__bowtie_phix:
    """Run bowtie phix"""
    input:
        [
            DECON_BOWTIE_PHIX / "bam" / f"{samples}.{library}_PhiX.bam"
            for samples, library in SAMPLES_LIBRARY
        ],
