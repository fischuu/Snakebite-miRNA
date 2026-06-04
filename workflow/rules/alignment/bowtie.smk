rule alignment__bowtie_mature_run:
    """
    Map non-tRNA and non-phiX samples against the mature.fa reference (bowtie).
    """
    input:
        reads = DECON_BOWTIE_PHIX / "unmapped" / "{samples}.{library}_PhiX_unmapped.fastq",
        index = REF_MIRBASE / "mature_basesAdjusted.fa.1.ebwt"
    output:
        mappedReads = ALIGN_BOWTIE_MATURE / "mapped" / "{samples}.{library}_mature_mapped.fastq" ,
        unmappedReads = ALIGN_BOWTIE_MATURE / "unmapped" / "{samples}.{library}_mature_unmapped.fastq" ,
        file = ALIGN_BOWTIE_MATURE / "bam" / "{samples}.{library}_mature.bam" ,
        wclmapped = ALIGN_BOWTIE_MATURE / "mapped" / "{samples}.{library}_mature_mapped.wcl" ,
        wclunmapped = ALIGN_BOWTIE_MATURE / "unmapped" / "{samples}.{library}_mature_unmapped.wcl" ,
        wccmapped = ALIGN_BOWTIE_MATURE / "mapped" / "{samples}.{library}_mature_mapped.wcc" ,
        wccunmapped = ALIGN_BOWTIE_MATURE / "unmapped" / "{samples}.{library}_mature_unmapped.wcc" 
    params:
        bamFolder=ALIGN_BOWTIE_MATURE / "bam",
        fastqFolder=ALIGN_BOWTIE_MATURE,
        index=REF_MIRBASE / "mature_basesAdjusted.fa",
        m=params["alignment"]["bowtie"]["m"],
        k=params["alignment"]["bowtie"]["k"]
    log:
        ALIGN_BOWTIE_MATURE / "logs" / "bowtie_map_mature.{samples}.{library}.log"
    benchmark:
        ALIGN_BOWTIE_MATURE / "benchmark" / "bowtie_map_mature.{samples}.{library}.benchmark.tsv"
    threads: esc("cpus", "alignment__bowtie_mature_run")
    resources:
        runtime=esc("runtime", "alignment__bowtie_mature_run"),
        mem_mb=esc("mem_mb", "alignment__bowtie_mature_run"),
        cpus_per_task=esc("cpus", "alignment__bowtie_mature_run"),
        slurm_partition=esc("partition", "alignment__bowtie_mature_run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'alignment__bowtie_mature_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("alignment__bowtie_mature_run"))
    container: docker["bowtie"]
    shell:"""
        mkdir -p {params.bamFolder}
        mkdir -p {params.fastqFolder}

        bowtie --best --strata --threads {threads} -k {params.k} -a -e 99999 --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input.reads}  | samtools view -bS | samtools sort - > {output.file} 2> {log};

        # I need to do this in case that there is a 0.00% alignment rate. In that case no files are created...
        touch {output.mappedReads}
        touch {output.unmappedReads}
        touch {output.file}
        
        wc -l {output.mappedReads} > {output.wclmapped}
        wc -l {output.unmappedReads} > {output.wclunmapped}
        sed -n '2~4p' {output.mappedReads}| wc -c > {output.wccmapped}
        sed -n '2~4p' {output.unmappedReads}| wc -c  > {output.wccunmapped}
    """

rule alignment__bowtie_hairpin_run:
    """
    Map non-tRNA and non-phiX samples against the hairpin.fa reference (bowtie).
    """
    input:
        reads = DECON_BOWTIE_PHIX / "unmapped" / "{samples}.{library}_PhiX_unmapped.fastq",
        index = REF_MIRBASE / "hairpin_basesAdjusted.fa.1.ebwt"
    output:
        mappedReads = ALIGN_BOWTIE_HAIRPIN / "mapped" / "{samples}.{library}_hairpin_mapped.fastq" ,
        unmappedReads = ALIGN_BOWTIE_HAIRPIN / "unmapped" / "{samples}.{library}_hairpin_unmapped.fastq" ,
        file = ALIGN_BOWTIE_HAIRPIN / "bam" / "{samples}.{library}_hairpin.bam" ,
        wclmapped = ALIGN_BOWTIE_HAIRPIN / "mapped" / "{samples}.{library}_hairpin_mapped.wcl" ,
        wclunmapped = ALIGN_BOWTIE_HAIRPIN / "unmapped" / "{samples}.{library}_hairpin_unmapped.wcl" ,
        wccmapped = ALIGN_BOWTIE_HAIRPIN / "mapped" / "{samples}.{library}_hairpin_mapped.wcc" ,
        wccunmapped = ALIGN_BOWTIE_HAIRPIN / "unmapped" / "{samples}.{library}_hairpin_unmapped.wcc" 
    params:
        bamFolder=ALIGN_BOWTIE_HAIRPIN / "bam",
        fastqFolder=ALIGN_BOWTIE_HAIRPIN,
        index=REF_MIRBASE / "hairpin_basesAdjusted.fa",
        m=params["alignment"]["bowtie"]["m"],
        k=params["alignment"]["bowtie"]["k"]
    log:
        ALIGN_BOWTIE_HAIRPIN / "logs" / "bowtie_map_hairpin.{samples}.{library}.log"
    benchmark:
        ALIGN_BOWTIE_HAIRPIN / "benchmark" / "bowtie_map_hairpin.{samples}.{library}.benchmark.tsv"
    threads: esc("cpus", "alignment__bowtie_hairpin_run")
    resources:
        runtime=esc("runtime", "alignment__bowtie_hairpin_run"),
        mem_mb=esc("mem_mb", "alignment__bowtie_hairpin_run"),
        cpus_per_task=esc("cpus", "alignment__bowtie_hairpin_run"),
        slurm_partition=esc("partition", "alignment__bowtie_hairpin_run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'alignment__bowtie_hairpin_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("alignment__bowtie_hairpin_run"))
    container: docker["bowtie"]
    shell:"""
        mkdir -p {params.bamFolder}
        mkdir -p {params.fastqFolder}

        bowtie --best --strata --threads {threads} -k {params.k} -a -e 99999 --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input.reads}  | samtools view -bS | samtools sort - > {output.file} 2> {log};

        # I need to do this in case that there is a 0.00% alignment rate. In that case no files are created...
        touch {output.mappedReads}
        touch {output.unmappedReads}
        touch {output.file}
        
        wc -l {output.mappedReads} > {output.wclmapped}
        wc -l {output.unmappedReads} > {output.wclunmapped}
        sed -n '2~4p' {output.mappedReads}| wc -c > {output.wccmapped}
        sed -n '2~4p' {output.unmappedReads}| wc -c  > {output.wccunmapped}
    """

rule alignment__bowtie_mature_species_run:
    """
    Map non-tRNA and non-phiX samples against the mature_species.fa reference (bowtie).
    """
    input:
        reads = DECON_BOWTIE_PHIX / "unmapped" / "{samples}.{library}_PhiX_unmapped.fastq",
        index = REF_MIRBASE / "mature_basesAdjusted_species.fa.1.ebwt"
    output:
        mappedReads = ALIGN_BOWTIE_MATURE_SPECIES / "mapped" / "{samples}.{library}_mature_species_mapped.fastq" ,
        unmappedReads = ALIGN_BOWTIE_MATURE_SPECIES / "unmapped" / "{samples}.{library}_mature_species_unmapped.fastq" ,
        file = ALIGN_BOWTIE_MATURE_SPECIES / "bam" / "{samples}.{library}_mature_species.bam" ,
        wclmapped = ALIGN_BOWTIE_MATURE_SPECIES / "mapped" / "{samples}.{library}_mature_species_mapped.wcl" ,
        wclunmapped = ALIGN_BOWTIE_MATURE_SPECIES / "unmapped" / "{samples}.{library}_mature_species_unmapped.wcl" ,
        wccmapped = ALIGN_BOWTIE_MATURE_SPECIES / "mapped" / "{samples}.{library}_mature_species_mapped.wcc" ,
        wccunmapped = ALIGN_BOWTIE_MATURE_SPECIES / "unmapped" / "{samples}.{library}_mature_species_unmapped.wcc" 
    params:
        bamFolder=ALIGN_BOWTIE_MATURE_SPECIES / "bam",
        fastqFolder=ALIGN_BOWTIE_MATURE_SPECIES,
        index=REF_MIRBASE / "mature_basesAdjusted_species.fa",
        m=params["alignment"]["bowtie"]["m"],
        k=params["alignment"]["bowtie"]["k"]
    log:
        ALIGN_BOWTIE_MATURE_SPECIES / "logs" / "bowtie_map_mature_species.{samples}.{library}.log"
    benchmark:
        ALIGN_BOWTIE_MATURE_SPECIES / "benchmark" / "bowtie_map_mature_species.{samples}.{library}.benchmark.tsv"
    threads: esc("cpus", "alignment__bowtie_mature_species_run")
    resources:
        runtime=esc("runtime", "alignment__bowtie_mature_species_run"),
        mem_mb=esc("mem_mb", "alignment__bowtie_mature_species_run"),
        cpus_per_task=esc("cpus", "alignment__bowtie_mature_species_run"),
        slurm_partition=esc("partition", "alignment__bowtie_mature_species_run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'alignment__bowtie_mature_species_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("alignment__bowtie_mature_species_run"))
    container: docker["bowtie"]
    shell:"""
        mkdir -p {params.bamFolder}
        mkdir -p {params.fastqFolder}

        bowtie --best --strata --threads {threads} -k {params.k} -a -e 99999 --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input.reads}  | samtools view -bS | samtools sort - > {output.file} 2> {log};

        # I need to do this in case that there is a 0.00% alignment rate. In that case no files are created...
        touch {output.mappedReads}
        touch {output.unmappedReads}
        touch {output.file}
        
        wc -l {output.mappedReads} > {output.wclmapped}
        wc -l {output.unmappedReads} > {output.wclunmapped}
        sed -n '2~4p' {output.mappedReads}| wc -c > {output.wccmapped}
        sed -n '2~4p' {output.unmappedReads}| wc -c  > {output.wccunmapped}
    """

rule alignment__bowtie_hairpin_species_run:
    """
    Map non-tRNA and non-phiX samples against the hairpin_species.fa reference (bowtie).
    """
    input:
        reads = DECON_BOWTIE_PHIX / "unmapped" / "{samples}.{library}_PhiX_unmapped.fastq",
        index = REF_MIRBASE / "hairpin_basesAdjusted_species.fa.1.ebwt"
    output:
        mappedReads = ALIGN_BOWTIE_HAIRPIN_SPECIES / "mapped" / "{samples}.{library}_hairpin_species_mapped.fastq" ,
        unmappedReads = ALIGN_BOWTIE_HAIRPIN_SPECIES / "unmapped" / "{samples}.{library}_hairpin_species_unmapped.fastq" ,
        file = ALIGN_BOWTIE_HAIRPIN_SPECIES / "bam" / "{samples}.{library}_hairpin_species.bam" ,
        wclmapped = ALIGN_BOWTIE_HAIRPIN_SPECIES / "mapped" / "{samples}.{library}_hairpin_species_mapped.wcl" ,
        wclunmapped = ALIGN_BOWTIE_HAIRPIN_SPECIES / "unmapped" / "{samples}.{library}_hairpin_species_unmapped.wcl" ,
        wccmapped = ALIGN_BOWTIE_HAIRPIN_SPECIES / "mapped" / "{samples}.{library}_hairpin_species_mapped.wcc" ,
        wccunmapped = ALIGN_BOWTIE_HAIRPIN_SPECIES / "unmapped" / "{samples}.{library}_hairpin_species_unmapped.wcc" 
    params:
        bamFolder=ALIGN_BOWTIE_HAIRPIN_SPECIES / "bam",
        fastqFolder=ALIGN_BOWTIE_HAIRPIN_SPECIES,
        index=REF_MIRBASE / "hairpin_basesAdjusted_species.fa",
        m=params["alignment"]["bowtie"]["m"],
        k=params["alignment"]["bowtie"]["k"]
    log:
        ALIGN_BOWTIE_HAIRPIN_SPECIES / "logs" / "bowtie_map_hairpin_species.{samples}.{library}.log"
    benchmark:
        ALIGN_BOWTIE_HAIRPIN_SPECIES / "benchmark" / "bowtie_map_hairpin_species.{samples}.{library}.benchmark.tsv"
    threads: esc("cpus", "alignment__bowtie_hairpin_species_run")
    resources:
        runtime=esc("runtime", "alignment__bowtie_hairpin_species_run"),
        mem_mb=esc("mem_mb", "alignment__bowtie_hairpin_species_run"),
        cpus_per_task=esc("cpus", "alignment__bowtie_hairpin_species_run"),
        slurm_partition=esc("partition", "alignment__bowtie_hairpin_species_run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'alignment__bowtie_hairpin_species_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("alignment__bowtie_hairpin_species_run"))
    container: docker["bowtie"]
    shell:"""
        mkdir -p {params.bamFolder}
        mkdir -p {params.fastqFolder}

        bowtie --best --strata --threads {threads} -k {params.k} -a -e 99999 --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input.reads}  | samtools view -bS | samtools sort - > {output.file} 2> {log};

        # I need to do this in case that there is a 0.00% alignment rate. In that case no files are created...
        touch {output.mappedReads}
        touch {output.unmappedReads}
        touch {output.file}
        
        wc -l {output.mappedReads} > {output.wclmapped}
        wc -l {output.unmappedReads} > {output.wclunmapped}
        sed -n '2~4p' {output.mappedReads}| wc -c > {output.wccmapped}
        sed -n '2~4p' {output.unmappedReads}| wc -c  > {output.wccunmapped}
    """
    
rule alignment__bowtie_reference_run:
    """
    Map non-tRNA, non-phiX samples against the reference genome (bowtie).
    """
    input:
        reads = DECON_BOWTIE_PHIX / "unmapped" / "{samples}.{library}_PhiX_unmapped.fastq",
        index = features["references"]["reference"] + ".1.ebwt"
    output:
        mappedReads = ALIGN_BOWTIE_REF / "mapped" / "{samples}.{library}_reference_mapped.fastq" ,
        unmappedReads = ALIGN_BOWTIE_REF / "unmapped" / "{samples}.{library}_reference_unmapped.fastq" ,
        file = ALIGN_BOWTIE_REF / "bam" / "{samples}.{library}_reference.bam" ,
        wclmapped = ALIGN_BOWTIE_REF / "mapped" / "{samples}.{library}_reference_mapped.wcl" ,
        wclunmapped = ALIGN_BOWTIE_REF / "unmapped" / "{samples}.{library}_reference_unmapped.wcl" ,
        wccmapped = ALIGN_BOWTIE_REF / "mapped" / "{samples}.{library}_reference_mapped.wcc" ,
        wccunmapped = ALIGN_BOWTIE_REF / "unmapped" / "{samples}.{library}_reference_unmapped.wcc" 
    params:
        bamFolder=ALIGN_BOWTIE_REF / "bam",
        fastqFolder=ALIGN_BOWTIE_REF,
        index=features["references"]["reference"],
        m=params["alignment"]["bowtie"]["m"],
        k=params["alignment"]["bowtie"]["k"]
    log:
        ALIGN_BOWTIE_REF / "logs" / "bowtie_map_reference.{samples}.{library}.log"
    benchmark:
        ALIGN_BOWTIE_REF / "benchmark" / "bowtie_map_reference.{samples}.{library}.benchmark.tsv"
    threads: esc("cpus", "alignment__bowtie_reference_run")
    resources:
        runtime=esc("runtime", "alignment__bowtie_reference_run"),
        mem_mb=esc("mem_mb", "alignment__bowtie_reference_run"),
        cpus_per_task=esc("cpus", "alignment__bowtie_reference_run"),
        slurm_partition=esc("partition", "alignment__bowtie_reference_run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'alignment__bowtie_reference_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("alignment__bowtie_reference_run"))
    container: docker["bowtie"]
    shell:"""
        mkdir -p {params.bamFolder}
        mkdir -p {params.fastqFolder}

        bowtie --best --strata --threads {threads} -k 50 -a -e 99999 --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input.reads} | samtools view -bS | samtools sort - > {output.file} 2> {log};

        # I need to do this in case that there is a 0.00% alignment rate. In that case no files are created...
        touch {output.mappedReads}
        touch {output.unmappedReads}
        touch {output.file}
        
        wc -l {output.mappedReads} > {output.wclmapped}
        wc -l {output.unmappedReads} > {output.wclunmapped}
        sed -n '2~4p' {output.mappedReads}| wc -c > {output.wccmapped}
        sed -n '2~4p' {output.unmappedReads}| wc -c  > {output.wccunmapped}
    """
    
if features["references"]["otherref"] == "":
    pass
else:
    rule alignment__bowtie_other_run:
        """
        Map non-tRNA and non-phiX samples against the other reference (bowtie).
        """
        input:
            reads = DECON_BOWTIE_PHIX / "unmapped" / "{samples}.{library}_PhiX_unmapped.fastq",
            index = features["references"]["otherref"] + ".1.ebwt"
        output:
            mappedReads = ALIGN_BOWTIE_OTHER / "mapped" / "{samples}.{library}_other_mapped.fastq" ,
            unmappedReads = ALIGN_BOWTIE_OTHER / "unmapped" / "{samples}.{library}_other_unmapped.fastq" ,
            file = ALIGN_BOWTIE_OTHER / "bam" / "{samples}.{library}_other.bam" ,
            wclmapped = ALIGN_BOWTIE_OTHER / "mapped" / "{samples}.{library}_other_mapped.wcl" ,
            wclunmapped = ALIGN_BOWTIE_OTHER / "unmapped" / "{samples}.{library}_other_unmapped.wcl" ,
            wccmapped = ALIGN_BOWTIE_OTHER / "mapped" / "{samples}.{library}_other_mapped.wcc" ,
            wccunmapped = ALIGN_BOWTIE_OTHER / "unmapped" / "{samples}.{library}_other_unmapped.wcc" 
        params:
            bamFolder=ALIGN_BOWTIE_OTHER / "bam",
            fastqFolder=ALIGN_BOWTIE_OTHER,
            index=features["references"]["otherref"],
            m=params["alignment"]["bowtie"]["m"],
            k=params["alignment"]["bowtie"]["k"]
        log:
            ALIGN_BOWTIE_OTHER / "logs" / "bowtie_map_other.{samples}.{library}.log"
        benchmark:
            ALIGN_BOWTIE_OTHER / "benchmark" / "bowtie_map_other.{samples}.{library}.benchmark.tsv"
        threads: esc("cpus", "alignment__bowtie_other_run")
        resources:
            runtime=esc("runtime", "alignment__bowtie_other_run"),
            mem_mb=esc("mem_mb", "alignment__bowtie_other_run"),
            cpus_per_task=esc("cpus", "alignment__bowtie_other_run"),
            slurm_partition=esc("partition", "alignment__bowtie_other_run"),
            gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'alignment__bowtie_other_run')['nvme']}",
            attempt=get_attempt,
        retries: len(get_escalation_order("alignment__bowtie_other_run"))
        container: docker["bowtie"]
        shell:"""
            mkdir -p {params.bamFolder}
            mkdir -p {params.fastqFolder}
    
            bowtie --best --strata --threads {threads} -k 50 -a -e 99999 --sam --al {output.mappedReads} --un {output.unmappedReads} {params.index} {input.reads}  | samtools view -bS | samtools sort - > {output.file} 2> {log};
    
            # I need to do this in case that there is a 0.00% alignment rate. In that case no files are created...
            touch {output.mappedReads}
            touch {output.unmappedReads}
            touch {output.file}
            
            wc -l {output.mappedReads} > {output.wclmapped}
            wc -l {output.unmappedReads} > {output.wclunmapped}
            sed -n '2~4p' {output.mappedReads}| wc -c > {output.wccmapped}
            sed -n '2~4p' {output.unmappedReads}| wc -c  > {output.wccunmapped}
        """
    

rule alignment__bowtie_mature:
    """Run bowtie mature"""
    input:
        [
            ALIGN_BOWTIE_MATURE / "bam" / f"{samples}.{library}_mature.bam"
            for samples, library in SAMPLES_LIBRARY
        ],

rule alignment__bowtie_hairpin:
    """Run bowtie hairpin"""
    input:
        [
            ALIGN_BOWTIE_HAIRPIN / "bam" / f"{samples}.{library}_hairpin.bam"
            for samples, library in SAMPLES_LIBRARY
        ],

rule alignment__bowtie_mature_species:
    """Run bowtie mature_species"""
    input:
        [
            ALIGN_BOWTIE_MATURE_SPECIES / "bam" / f"{samples}.{library}_mature_species.bam"
            for samples, library in SAMPLES_LIBRARY
        ],

rule alignment__bowtie_hairpin_species:
    """Run bowtie hairpin_species"""
    input:
        [
            ALIGN_BOWTIE_HAIRPIN_SPECIES / "bam" / f"{samples}.{library}_hairpin_species.bam"
            for samples, library in SAMPLES_LIBRARY
        ],
        
rule alignment__bowtie_reference:
    """Run bowtie reference"""
    input:
        [
            ALIGN_BOWTIE_REF / "bam" / f"{samples}.{library}_reference.bam"
            for samples, library in SAMPLES_LIBRARY
        ],
        
if features["references"]["otherref"] == "":
    pass
else:
    rule alignment__bowtie_other:
        """Run bowtie reference"""
        input:
            [
                ALIGN_BOWTIE_OTHER / "bam" / f"{samples}.{library}_other.bam"
                for samples, library in SAMPLES_LIBRARY
            ],
