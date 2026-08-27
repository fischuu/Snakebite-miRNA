rule align__star__mature:
    """Map PhiX-unmapped reads against mature.fa (STAR)"""
    input:
        index=rules.reference__star_index__mature.output,
        fastq=rules.decontaminate__bowtie__phix.output.unmapped,
    output:
        bam=MATURE_STAR / "{sample_id}_mature_star.bam",
        logdir=directory(MATURE_STAR / "{sample_id}"),
    log:
        MATURE_STAR / "star.{sample_id}.log"
    benchmark:
        MATURE_STAR / "benchmark/star.{sample_id}.tsv"
    params:
        index=str(STAR_INDEX_MATURE),
        tmpdir=str(STAR_ALIGN_TMP / "mature" / "{sample_id}"),
        limit_bam_sort_ram=params["align"]["star"]["limit_bam_sort_ram"],
    threads: esc("cpus", "align__star__mature")
    resources:
        runtime=esc("runtime", "align__star__mature"),
        mem_mb=esc("mem_mb", "align__star__mature"),
        cpus_per_task=esc("cpus", "align__star__mature"),
        slurm_partition=esc("partition", "align__star__mature"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__star__mature')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__star__mature"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {output.logdir} {params.tmpdir}
        cd {params.tmpdir}/..
        STAR --genomeDir {params.index} --outTmpDir {params.tmpdir} \
             --readFilesIn {input.fastq} \
             --outFilterMismatchNoverLmax 0.01 --outFilterMatchNmin 16 \
             --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
             --outFilterMultimapNmax 50 --chimSegmentMin 20 --chimOutType WithinBAM \
             --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
             --outMultimapperOrder Random --runThreadN {threads} \
             --limitBAMsortRAM {params.limit_bam_sort_ram} \
             --outFileNamePrefix {wildcards.sample_id}_
        mv {wildcards.sample_id}_Aligned.sortedByCoord.out.bam {output.bam}
        mv {wildcards.sample_id}_Log.final.out {wildcards.sample_id}_Log.progress.out \
           {wildcards.sample_id}_Log.out {wildcards.sample_id}_SJ.out.tab {output.logdir}
        """

rule align__samtools__star_mature_flagstat:
    """Report mapping stats and extract mapped/unmapped FASTQ from the mature STAR alignment (samtools)"""
    input:
        rules.align__star__mature.output.bam
    output:
        flagstat=MATURE_STAR / "{sample_id}_mature.flagstat",
        stats=MATURE_STAR / "{sample_id}_mature.stats",
        mapped=MATURE_STAR / "mapped" / "{sample_id}_mature_mapped_star.fastq",
        unmapped=MATURE_STAR / "unmapped" / "{sample_id}_mature_unmapped_star.fastq",
    log:
        MATURE_STAR / "samtools_flagstat.{sample_id}.log"
    benchmark:
        MATURE_STAR / "benchmark/samtools_flagstat.{sample_id}.tsv"
    threads: esc("cpus", "align__samtools__star_mature_flagstat")
    resources:
        runtime=esc("runtime", "align__samtools__star_mature_flagstat"),
        mem_mb=esc("mem_mb", "align__samtools__star_mature_flagstat"),
        cpus_per_task=esc("cpus", "align__samtools__star_mature_flagstat"),
        slurm_partition=esc("partition", "align__samtools__star_mature_flagstat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__samtools__star_mature_flagstat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__samtools__star_mature_flagstat"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats}
        samtools bam2fq -F 4 {input} > {output.mapped}
        samtools bam2fq -f 4 {input} > {output.unmapped}
        touch {output.mapped} {output.unmapped}
        """

rule align__star__mature_species:
    """Map PhiX-unmapped reads against the species-filtered mature.fa (STAR)"""
    input:
        index=rules.reference__star_index__mature_species.output,
        fastq=rules.decontaminate__bowtie__phix.output.unmapped,
    output:
        bam=MATURE_SPECIES_STAR / "{sample_id}_mature_species_star.bam",
        logdir=directory(MATURE_SPECIES_STAR / "{sample_id}"),
    log:
        MATURE_SPECIES_STAR / "star.{sample_id}.log"
    benchmark:
        MATURE_SPECIES_STAR / "benchmark/star.{sample_id}.tsv"
    params:
        index=str(STAR_INDEX_MATURE_SPECIES),
        tmpdir=str(STAR_ALIGN_TMP / "mature_species" / "{sample_id}"),
        limit_bam_sort_ram=params["align"]["star"]["limit_bam_sort_ram"],
    threads: esc("cpus", "align__star__mature_species")
    resources:
        runtime=esc("runtime", "align__star__mature_species"),
        mem_mb=esc("mem_mb", "align__star__mature_species"),
        cpus_per_task=esc("cpus", "align__star__mature_species"),
        slurm_partition=esc("partition", "align__star__mature_species"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__star__mature_species')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__star__mature_species"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {output.logdir} {params.tmpdir}
        cd {params.tmpdir}/..
        STAR --genomeDir {params.index} --outTmpDir {params.tmpdir} \
             --readFilesIn {input.fastq} \
             --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 \
             --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
             --outFilterMultimapNmax 50 --chimSegmentMin 20 --chimOutType WithinBAM \
             --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
             --outMultimapperOrder Random --runThreadN {threads} \
             --limitBAMsortRAM {params.limit_bam_sort_ram} \
             --outFileNamePrefix {wildcards.sample_id}_
        mv {wildcards.sample_id}_Aligned.sortedByCoord.out.bam {output.bam}
        mv {wildcards.sample_id}_Log.final.out {wildcards.sample_id}_Log.progress.out \
           {wildcards.sample_id}_Log.out {wildcards.sample_id}_SJ.out.tab {output.logdir}
        """

rule align__samtools__star_mature_species_flagstat:
    """Report mapping stats and extract mapped/unmapped FASTQ from the mature_species STAR alignment"""
    input:
        rules.align__star__mature_species.output.bam
    output:
        flagstat=MATURE_SPECIES_STAR / "{sample_id}_mature_species.flagstat",
        stats=MATURE_SPECIES_STAR / "{sample_id}_mature_species.stats",
        mapped=MATURE_SPECIES_STAR / "mapped" / "{sample_id}_mature_species_mapped_star.fastq",
        unmapped=MATURE_SPECIES_STAR / "unmapped" / "{sample_id}_mature_species_unmapped_star.fastq",
    log:
        MATURE_SPECIES_STAR / "samtools_flagstat.{sample_id}.log"
    benchmark:
        MATURE_SPECIES_STAR / "benchmark/samtools_flagstat.{sample_id}.tsv"
    threads: esc("cpus", "align__samtools__star_mature_species_flagstat")
    resources:
        runtime=esc("runtime", "align__samtools__star_mature_species_flagstat"),
        mem_mb=esc("mem_mb", "align__samtools__star_mature_species_flagstat"),
        cpus_per_task=esc("cpus", "align__samtools__star_mature_species_flagstat"),
        slurm_partition=esc("partition", "align__samtools__star_mature_species_flagstat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__samtools__star_mature_species_flagstat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__samtools__star_mature_species_flagstat"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats}
        samtools bam2fq -F 4 {input} > {output.mapped}
        samtools bam2fq -f 4 {input} > {output.unmapped}
        touch {output.mapped} {output.unmapped}
        """

rule align__star__hairpin:
    """Map PhiX-unmapped reads against hairpin.fa (STAR)"""
    input:
        index=rules.reference__star_index__hairpin.output,
        fastq=rules.decontaminate__bowtie__phix.output.unmapped,
    output:
        bam=HAIRPIN_STAR / "{sample_id}_hairpin_star.bam",
        logdir=directory(HAIRPIN_STAR / "{sample_id}"),
    log:
        HAIRPIN_STAR / "star.{sample_id}.log"
    benchmark:
        HAIRPIN_STAR / "benchmark/star.{sample_id}.tsv"
    params:
        index=str(STAR_INDEX_HAIRPIN),
        tmpdir=str(STAR_ALIGN_TMP / "hairpin" / "{sample_id}"),
        limit_bam_sort_ram=params["align"]["star"]["limit_bam_sort_ram"],
    threads: esc("cpus", "align__star__hairpin")
    resources:
        runtime=esc("runtime", "align__star__hairpin"),
        mem_mb=esc("mem_mb", "align__star__hairpin"),
        cpus_per_task=esc("cpus", "align__star__hairpin"),
        slurm_partition=esc("partition", "align__star__hairpin"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__star__hairpin')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__star__hairpin"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {output.logdir} {params.tmpdir}
        cd {params.tmpdir}/..
        STAR --genomeDir {params.index} --outTmpDir {params.tmpdir} \
             --readFilesIn {input.fastq} \
             --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 \
             --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
             --outFilterMultimapNmax 50 --chimSegmentMin 20 --chimOutType WithinBAM \
             --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
             --outMultimapperOrder Random --runThreadN {threads} \
             --limitBAMsortRAM {params.limit_bam_sort_ram} \
             --outFileNamePrefix {wildcards.sample_id}_
        mv {wildcards.sample_id}_Aligned.sortedByCoord.out.bam {output.bam}
        mv {wildcards.sample_id}_Log.final.out {wildcards.sample_id}_Log.progress.out \
           {wildcards.sample_id}_Log.out {wildcards.sample_id}_SJ.out.tab {output.logdir}
        """

rule align__samtools__star_hairpin_flagstat:
    """Report mapping stats and extract mapped/unmapped FASTQ from the hairpin STAR alignment"""
    input:
        rules.align__star__hairpin.output.bam
    output:
        flagstat=HAIRPIN_STAR / "{sample_id}_hairpin.flagstat",
        stats=HAIRPIN_STAR / "{sample_id}_hairpin.stats",
        mapped=HAIRPIN_STAR / "mapped" / "{sample_id}_hairpin_mapped_star.fastq",
        unmapped=HAIRPIN_STAR / "unmapped" / "{sample_id}_hairpin_unmapped_star.fastq",
    log:
        HAIRPIN_STAR / "samtools_flagstat.{sample_id}.log"
    benchmark:
        HAIRPIN_STAR / "benchmark/samtools_flagstat.{sample_id}.tsv"
    threads: esc("cpus", "align__samtools__star_hairpin_flagstat")
    resources:
        runtime=esc("runtime", "align__samtools__star_hairpin_flagstat"),
        mem_mb=esc("mem_mb", "align__samtools__star_hairpin_flagstat"),
        cpus_per_task=esc("cpus", "align__samtools__star_hairpin_flagstat"),
        slurm_partition=esc("partition", "align__samtools__star_hairpin_flagstat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__samtools__star_hairpin_flagstat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__samtools__star_hairpin_flagstat"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats}
        samtools bam2fq -F 4 {input} > {output.mapped}
        samtools bam2fq -f 4 {input} > {output.unmapped}
        touch {output.mapped} {output.unmapped}
        """

rule align__star__hairpin_species:
    """Map mature_species STAR-unmapped reads against the species-filtered hairpin.fa (STAR)"""
    input:
        index=rules.reference__star_index__hairpin_species.output,
        fastq=rules.align__samtools__star_mature_species_flagstat.output.unmapped,
    output:
        bam=HAIRPIN_SPECIES_STAR / "{sample_id}_hairpin_species_star.bam",
        logdir=directory(HAIRPIN_SPECIES_STAR / "{sample_id}"),
    log:
        HAIRPIN_SPECIES_STAR / "star.{sample_id}.log"
    benchmark:
        HAIRPIN_SPECIES_STAR / "benchmark/star.{sample_id}.tsv"
    params:
        index=str(STAR_INDEX_HAIRPIN_SPECIES),
        tmpdir=str(STAR_ALIGN_TMP / "hairpin_species" / "{sample_id}"),
        limit_bam_sort_ram=params["align"]["star"]["limit_bam_sort_ram"],
    threads: esc("cpus", "align__star__hairpin_species")
    resources:
        runtime=esc("runtime", "align__star__hairpin_species"),
        mem_mb=esc("mem_mb", "align__star__hairpin_species"),
        cpus_per_task=esc("cpus", "align__star__hairpin_species"),
        slurm_partition=esc("partition", "align__star__hairpin_species"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__star__hairpin_species')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__star__hairpin_species"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {output.logdir} {params.tmpdir}
        cd {params.tmpdir}/..
        STAR --genomeDir {params.index} --outTmpDir {params.tmpdir} \
             --readFilesIn {input.fastq} \
             --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 \
             --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
             --outFilterMultimapNmax 50 --chimSegmentMin 20 --chimOutType WithinBAM \
             --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
             --outMultimapperOrder Random --runThreadN {threads} \
             --limitBAMsortRAM {params.limit_bam_sort_ram} \
             --outFileNamePrefix {wildcards.sample_id}_
        mv {wildcards.sample_id}_Aligned.sortedByCoord.out.bam {output.bam}
        mv {wildcards.sample_id}_Log.final.out {wildcards.sample_id}_Log.progress.out \
           {wildcards.sample_id}_Log.out {wildcards.sample_id}_SJ.out.tab {output.logdir}
        """

rule align__samtools__star_hairpin_species_flagstat:
    """Report mapping stats and extract mapped/unmapped FASTQ from the hairpin_species STAR alignment"""
    input:
        rules.align__star__hairpin_species.output.bam
    output:
        flagstat=HAIRPIN_SPECIES_STAR / "{sample_id}_hairpin_species.flagstat",
        stats=HAIRPIN_SPECIES_STAR / "{sample_id}_hairpin_species.stats",
        mapped=HAIRPIN_SPECIES_STAR / "mapped" / "{sample_id}_hairpin_species_mapped_star.fastq",
        unmapped=HAIRPIN_SPECIES_STAR / "unmapped" / "{sample_id}_hairpin_species_unmapped_star.fastq",
    log:
        HAIRPIN_SPECIES_STAR / "samtools_flagstat.{sample_id}.log"
    benchmark:
        HAIRPIN_SPECIES_STAR / "benchmark/samtools_flagstat.{sample_id}.tsv"
    threads: esc("cpus", "align__samtools__star_hairpin_species_flagstat")
    resources:
        runtime=esc("runtime", "align__samtools__star_hairpin_species_flagstat"),
        mem_mb=esc("mem_mb", "align__samtools__star_hairpin_species_flagstat"),
        cpus_per_task=esc("cpus", "align__samtools__star_hairpin_species_flagstat"),
        slurm_partition=esc("partition", "align__samtools__star_hairpin_species_flagstat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__samtools__star_hairpin_species_flagstat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__samtools__star_hairpin_species_flagstat"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats}
        samtools bam2fq -F 4 {input} > {output.mapped}
        samtools bam2fq -f 4 {input} > {output.unmapped}
        touch {output.mapped} {output.unmapped}
        """

rule align__star__genome:
    """Map hairpin_species STAR-unmapped reads (non-miRBase) against the host reference genome (STAR)"""
    input:
        index=rules.reference__star_index__genome.output,
        fastq=rules.align__samtools__star_hairpin_species_flagstat.output.unmapped,
    output:
        bam=GENOME_STAR / "{sample_id}_reference_star.bam",
        logdir=directory(GENOME_STAR / "{sample_id}"),
    log:
        GENOME_STAR / "star.{sample_id}.log"
    benchmark:
        GENOME_STAR / "benchmark/star.{sample_id}.tsv"
    params:
        index=str(STAR_INDEX_GENOME),
        tmpdir=str(STAR_ALIGN_TMP / "genome" / "{sample_id}"),
        limit_bam_sort_ram=params["align"]["star"]["limit_bam_sort_ram"],
    threads: esc("cpus", "align__star__genome")
    resources:
        runtime=esc("runtime", "align__star__genome"),
        mem_mb=esc("mem_mb", "align__star__genome"),
        cpus_per_task=esc("cpus", "align__star__genome"),
        slurm_partition=esc("partition", "align__star__genome"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__star__genome')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__star__genome"))
    container:
        docker["star"]
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {output.logdir} {params.tmpdir}
        cd {params.tmpdir}/..
        STAR --genomeDir {params.index} --outTmpDir {params.tmpdir} \
             --readFilesIn {input.fastq} \
             --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 \
             --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
             --outFilterMultimapNmax 50 --chimSegmentMin 20 --chimOutType WithinBAM \
             --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
             --outMultimapperOrder Random --runThreadN {threads} \
             --limitBAMsortRAM {params.limit_bam_sort_ram} \
             --outFileNamePrefix {wildcards.sample_id}_
        mv {wildcards.sample_id}_Aligned.sortedByCoord.out.bam {output.bam}
        mv {wildcards.sample_id}_Log.final.out {wildcards.sample_id}_Log.progress.out \
           {wildcards.sample_id}_Log.out {wildcards.sample_id}_SJ.out.tab {output.logdir}
        """

rule align__samtools__remove_secondary:
    """Keep only primary alignments from the genome STAR bam (samtools)"""
    input:
        rules.align__star__genome.output.bam
    output:
        GENOME_STAR / "{sample_id}_reference_star_primary.bam"
    log:
        GENOME_STAR / "samtools_remove_secondary.{sample_id}.log"
    benchmark:
        GENOME_STAR / "benchmark/samtools_remove_secondary.{sample_id}.tsv"
    threads: esc("cpus", "align__samtools__remove_secondary")
    resources:
        runtime=esc("runtime", "align__samtools__remove_secondary"),
        mem_mb=esc("mem_mb", "align__samtools__remove_secondary"),
        cpus_per_task=esc("cpus", "align__samtools__remove_secondary"),
        slurm_partition=esc("partition", "align__samtools__remove_secondary"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__samtools__remove_secondary')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__samtools__remove_secondary"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools view -bq 1 {input} > {output}
        samtools index {output}
        """

rule align__samtools__star_genome_flagstat:
    """Report mapping stats and extract mapped/unmapped FASTQ from the genome STAR alignment"""
    input:
        rules.align__star__genome.output.bam
    output:
        flagstat=GENOME_STAR / "{sample_id}_reference.flagstat",
        stats=GENOME_STAR / "{sample_id}_reference.stats",
        mapped=GENOME_STAR / "mapped" / "{sample_id}_reference_mapped_star.fastq",
        unmapped=GENOME_STAR / "unmapped" / "{sample_id}_reference_unmapped_star.fastq",
    log:
        GENOME_STAR / "samtools_flagstat.{sample_id}.log"
    benchmark:
        GENOME_STAR / "benchmark/samtools_flagstat.{sample_id}.tsv"
    threads: esc("cpus", "align__samtools__star_genome_flagstat")
    resources:
        runtime=esc("runtime", "align__samtools__star_genome_flagstat"),
        mem_mb=esc("mem_mb", "align__samtools__star_genome_flagstat"),
        cpus_per_task=esc("cpus", "align__samtools__star_genome_flagstat"),
        slurm_partition=esc("partition", "align__samtools__star_genome_flagstat"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'align__samtools__star_genome_flagstat')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("align__samtools__star_genome_flagstat"))
    container:
        docker["samtools"]
    shell:
        """
        exec > {log} 2>&1
        samtools flagstat {input} > {output.flagstat}
        samtools stats {input} > {output.stats}
        samtools bam2fq -F 4 {input} > {output.mapped}
        samtools bam2fq -f 4 {input} > {output.unmapped}
        touch {output.mapped} {output.unmapped}
        """
