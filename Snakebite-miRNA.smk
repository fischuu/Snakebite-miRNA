import pandas as pd
from snakemake.utils import validate, min_version
from multiprocessing import cpu_count
import glob
import re
import os, sys
import yaml

##### Snakebite miRNA pipeline #####
##### Daniel Fischer (daniel.fischer@luke.fi)
##### Natural Resources Institute Finland (Luke)
##### Version: 0.6.2
version = "0.6.2"

##### set minimum snakemake version #####
min_version("6.0")

##### Sample sheets #####

##### load config and sample sheets #####

samplesheet = pd.read_table(config["samplesheet-file"]).set_index("rawsample", drop=False)
rawsamples=list(samplesheet.rawsample)
samples=list(set(list(samplesheet.sample_name)))
lane=list(samplesheet.lane)

wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples)
    
    ##### input function definitions ######

def get_lane(wildcards):
    output = samplesheet.loc[wildcards.rawsamples][["lane"]]
    return output.tolist()

def get_sample(wildcards):
    output = samplesheet.loc[wildcards.rawsamples][["sample_name"]]
    return output.tolist()

def get_raw_input_fastqs(wildcards):
    reads = samplesheet.loc[wildcards.rawsamples][["read1", "read2"]]
    path = config["rawdata-folder"]
    output = [path + x for x in reads]
    return output

def get_raw_input_read1(wildcards):
    reads = samplesheet.loc[wildcards.rawsamples][["read1"]]
    path = config["rawdata-folder"]
    output = [path + "/" + x for x in reads]
    return output

def get_fastq_for_concatenating_read1(wildcards):
    rs = samplesheet.loc[samplesheet["sample_name"] == wildcards.samples]["rawsample"]
    path = config["project-folder"] + "/FASTQ/TRIMMED/"
    output = [path + x for x in rs + "_trimmed.fastq.gz"]
    return output   

#### CONTINUE FROM HERE TO ADD PIPE CONFIG ONTO THE FILE
#if '--configfile' in sys.argv:
#    i = sys.argv.index('--configfile')
#    config["pipeline-config"] = sys.argv[i + 1]
#else:
#    config["pipeline-config"] = ""
#    
#if '--cluster-config' in sys.argv:
#    i = sys.argv.index('--cluster-config')
#    config["server-config"] = sys.argv[i + 1]
#else:
#    config["server-config"] = ""

##### Extract the cluster resource requests from the server config #####
cluster=dict()
if os.path.exists(config["server-config"]):
    with open(config["server-config"]) as yml:
        cluster = yaml.load(yml, Loader=yaml.FullLoader)

##### Complete the input configuration
config["tRNAIndex"] = config["tRNARef"]
config["phixIndex"] = config["phixRef"]
config["matureIndex"] =  "%s/References/mature_basesAdjusted.fa" % (config["project-folder"])
config["matureSpeciesIndex"] =  "%s/References/mature_basesAdjusted_species.fa" % (config["project-folder"])
config["hairpinIndex"] =  "%s/References/hairpin_basesAdjusted.fa" % (config["project-folder"])
config["referenceIndex"] = config["reference"]
config["matureSTARIndex"] = config["starbase"]+"/Mature"
config["hairpinSTARIndex"] = config["starbase"]+"/Hairpin"
config["matureSpeciesSTARIndex"] = config["starbase"]+"/MatureSpecies"
config["hairpinSpeciesSTARIndex"] = config["starbase"]+"/HairpinSpecies"
config["referenceSTARIndex"] = config["starbase"]+"/Reference"
config["report-script"] = config["pipeline-folder"]+"/scripts/workflow-report.Rmd"

##### Singularity container #####
config["singularity"] = {}
config["singularity"]["bedtools"] = "docker://fischuu/bedtools:2.30-0.1"
config["singularity"]["bowtie"] = "docker://fischuu/bowtie:1.2.2-0.3"
config["singularity"]["cutadapt"] = "docker://fischuu/cutadapt:2.8-0.3"
config["singularity"]["fastqc"] = "docker://fischuu/gbs:0.2"
config["singularity"]["multiqc"] = "docker://fischuu/gbs:0.2"
config["singularity"]["samtools"] = "docker://fischuu/samtools:1.9-0.2"
config["singularity"]["star"] = "docker://fischuu/star:2.7.3a-0.2"
config["singularity"]["subread"] = "docker://fischuu/subread:2.0.1-0.1"
config["singularity"]["seqkit"] = "docker://fischuu/seqkit:2.1.0-0.2"
config["singularity"]["reporting"] = "docker://fischuu/r-gbs:4.1.2-0.4"
config["singularity"]["stringtie"] = "docker://fischuu/stringtie:2.2.1-0.1"

##### Apply pre-configuration settings #####
if config["params"]["protocol"] == 'illumina':
    config["params"]["cutadapt"]["adapter"] = "TGGAATTCTCGGGTGCCAAGG"
    config["params"]["cutadapt"]["fiveprimetrim"] = 0 
    config["params"]["cutadapt"]["threeprimetrim"] = 0
elif config["params"]["protocol"] == 'nextflex':
    config["params"]["cutadapt"]["adapter"] = "TGGAATTCTCGGGTGCCAAGG"
    config["params"]["cutadapt"]["fiveprimetrim"] = 4
    config["params"]["cutadapt"]["threeprimetrim"] = 4
elif config["params"]["protocol"] == 'qiagen':
    config["params"]["cutadapt"]["adapter"] = "AACTGTAGGCACCATCAAT"
    config["params"]["cutadapt"]["fiveprimetrim"] = 0
    config["params"]["cutadapt"]["threeprimetrim"] = 0
elif config["params"]["protocol"] == 'lexogen':
    config["params"]["cutadapt"]["adapter"] = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
    config["params"]["cutadapt"]["fiveprimetrim"] = 0
    config["params"]["cutadapt"]["threeprimetrim"] = 0
    
##### Input constraints #####
wildcard_constraints:
    rawsamples="|".join(rawsamples),
    samples="|".join(samples)

##### Print the welcome screen #####
print("#################################################################################")
print("##### Welcome to the Snakemake miRNA pipeline")
print("##### Version: "+version)
print("#####")
print("##### Pipeline configuration")
print("##### --------------------------------")
print("##### project-folder  : "+config["project-folder"])
print("##### pipeline-folder : "+config["pipeline-folder"])
print("##### server-config   : "+config["server-config"])
print("##### pipeline-config : "+config["pipeline-config"])
print("##### species-id      : "+config["params"]["species-id"])
print("##### protocol        : "+config["params"]["protocol"])
print("#####")
print("##### Trimming configuration")
print("##### --------------------------------")
print("##### adapter         : "+ config["params"]["cutadapt"]["adapter"])
print("##### 5'-trim bases   : "+ str(config["params"]["cutadapt"]["fiveprimetrim"]))
print("##### 3'-trim bases   : "+ str(config["params"]["cutadapt"]["threeprimetrim"]))
print("#####")
print("##### Singularity configuration")
print("##### --------------------------------")
print("##### bedtools        : "+config["singularity"]["bedtools"])
print("##### bowtie          : "+config["singularity"]["bowtie"])
print("##### cutadapt        : "+config["singularity"]["cutadapt"])
print("##### fastqc          : "+config["singularity"]["fastqc"])
print("##### multiqc         : "+config["singularity"]["multiqc"])
print("##### samtools        : "+config["singularity"]["samtools"])
print("##### star            : "+config["singularity"]["star"])
print("##### subread         : "+config["singularity"]["subread"])
print("##### seqkit          : "+config["singularity"]["seqkit"])
print("##### stringtie       : "+config["singularity"]["stringtie"])
print("##### reporting       : "+config["singularity"]["reporting"])
print("#################################################################################")

##### run complete pipeline #####

rule all:
    input:
      # Module 1: QUALITY CHECKS, TRIMMING and PREPARATIONS
        "%s/QC/RAW/multiqc/" % (config["project-folder"]),
        "%s/QC/TRIMMED/multiqc/" % (config["project-folder"]),
        "%s/QC/CONCATENATED/multiqc/" % (config["project-folder"]),
        "%s/QC/tRNA/multiqc_unmapped/" % (config["project-folder"]),
        "%s/QC/tRNA/multiqc_mapped/" % (config["project-folder"]),
        "%s/QC/PhiX/multiqc_unmapped/" % (config["project-folder"]),
        "%s/QC/PhiX/multiqc_mapped/" % (config["project-folder"]),
        "%s/QC/Mature/multiqc_unmapped/" % (config["project-folder"]),
        "%s/QC/Mature/multiqc_mapped/" % (config["project-folder"]),
        "%s/QC/Hairpin/multiqc_unmapped/" % (config["project-folder"]),
        "%s/QC/Hairpin/multiqc_mapped/" % (config["project-folder"]),
        expand("%s/STATS/BOWTIE/tRNA/{samples}_tRNA.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/STATS/BOWTIE/PhiX/{samples}_PhiX.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/STATS/BOWTIE/Mature/{samples}_mature.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/STATS/BOWTIE/Hairpin/{samples}_hairpin.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/STATS/STAR/Mature/{samples}_mature.flagstat" % (config["project-folder"]), samples=samples),
        expand("%s/STATS/STAR/Reference/{samples}_reference.flagstat" % (config["project-folder"]), samples=samples),
      # ALIGNMENTS
        expand("%s/BAM/BOWTIE/tRNA/{samples}_tRNA.sorted.bam.bai" % (config["project-folder"]),samples=samples),
      # QUANTIFICATION
        expand("%s/FASTA/STAR/Reference_softclipped/{samples}_reference_softclipped.fasta.gz" % (config["project-folder"]), samples=samples),
        expand("%s/QUANTIFICATION/STAR/Reference/{samples}_star_reference_fc.txt" % (config["project-folder"]), samples=samples),
        expand("%s/QUANTIFICATION/BOWTIE/Mature/{samples}_bowtie_mature_seqkit.txt" % (config["project-folder"]), samples=samples),
        expand("%s/QUANTIFICATION/STAR/Reference/{samples}_star_reference_exon_fc.txt" % (config["project-folder"]), samples=samples),
      # NOVEL MIRNA
        expand("%s/QUANTIFICATION/STAR/Novel_genes/{samples}_star_novelMirna_bedtools.txt" % (config["project-folder"]), samples=samples),
      # REPORTING
        "%s/pipelineReport.html" % (config["project-folder"])
#      # ANALYSIS - Conservation
#        "%s/FASTA/ARS/ARS-UCD1.2.95.mirna.fa" % (config["project-folder"]),
#        "%s/FASTA/hg38/hg38.95.mirna.fa" % (config["project-folder"]),
#        "%s/Conservation/mirnaSeq.sam" % (config["project-folder"]),
#      # ANALYSIS - Alignments
#        expand("%s/BAM/STAR_softclipped/{samples}_star_softclipped.bam" % (config["project-folder"]), samples=samples),
#        expand("%s/BAM/STAR/{samples}_star.bam" % (config["project-folder"]), samples=samples),
#        expand("%s/FASTA/STAR/{samples}_softclipped.fasta.gz" % (config["project-folder"]), samples=samples),
#        expand("%s/BAM/ars/{samples}_ars_softclipped.bam" % (config["project-folder"]), samples=samples),
#        expand("%s/GTF/STAR/{samples}_star_fc.txt" % (config["project-folder"]), samples=samples),
#        expand("%s/GTF/STAR_softclipped/{samples}_softclipped_fc.txt" % (config["project-folder"]), samples=samples),
#        "%s/finalReport.html" % (config["project-folder"])
#        #expand("%s/GTF/STAR_cufflinks/{samples}_star.gtf" % (config["project-folder"]), samples=samples)
#        "%s/References/gtf_merged_star.gtf" % (config["project-folder"])
##      #  "%s/PRODIGAL/ars_unmapped/final.contigs.prodigal.gtf" % (config["project-folder"]),
##      #  "%s/kraken_ars_unmapped/ars_unmapped_taxonomy.report" % (config["project-folder"])
rule decontamination:
    input:
        expand("%s/FASTQ/tRNA/mapped/{samples}_tRNA_mapped.fastq" % (config["project-folder"]), samples=samples),
        expand("%s/FASTQ/PhiX/mapped/{samples}_PhiX_mapped.fastq" % (config["project-folder"]), samples=samples)
#### setup report #####
#
#report: "report/workflow.rst"
#
##### load rules #####

include: "rules/Module0-Preparations"
include: "rules/Module1-QC"
include: "rules/Module2-Preprocessing"
include: "rules/Module3-ContaminationRemoval"
include: "rules/Module4-Alignments"
include: "rules/Module5-Quantification"
include: "rules/Module6-Reporting"
include: "rules/Module7-NovelMirna"
#include: "rules/bedtools_extract_miRNA_fasta_ars.smk"
#include: "rules/bedtools_extract_miRNA_fasta_hg38.smk"
#include: "rules/bowtie_map_ars_mirna_vs_hg38.smk"
#include: "rules/cutadapt_trimming.smk"
#include: "rules/fastqc_quality_control_trimmedData.smk"
#include: "rules/multiqc_quality_control_trimmedData.smk"
#include: "rules/bash_catFastq.smk"
#include: "rules/fastqc_quality_control_concatenatedData.smk"
#include: "rules/multiqc_quality_control_concatenatedData.smk"
#include: "rules/bowtie_map_concatenated_reads_tRNA.smk"
#include: "rules/fastqc_quality_control_tRNA.smk"
#include: "rules/multiqc_quality_control_tRNA.smk"
#include: "rules/bowtie_map_nonTRNA_reads_phix.smk"
#include: "rules/bowtie_map_nonPhiX_reads_ars.smk"
#include: "rules/star_map_nonPhiX_reads_ars.smk"
#include: "rules/bash_get_softclipped_fastq.smk"
#include: "rules/bash_filter_softclipped_fastq.smk"
#include: "rules/bowtie_map_softclipped_reads_ars.smk"
#include: "rules/sam_to_sortedBam_ars_softclipped.smk"
#include: "rules/fastqc_quality_control_ars.smk"
#include: "rules/multiqc_quality_control_ars.smk"
#include: "rules/sam_to_sortedBam_ars.smk"
#include: "rules/cufflinks_quantify_star.smk"
#include: "rules/cuffmerge_compose_merge_quant_ars.smk"
#include: "rules/cuffmerge_merge_quant_ars.smk"
#include: "rules/cufflinks_quantify_merged_ars.smk"
#include: "rules/featureCounts_quantify_merged_ars.smk"
#include: "rules/bash_compose_assemble_ars_unmapped_reads.smk"
#include: "rules/fastx_trim_unmapped_ars.smk"
#include: "rules/bowtie_map_unmapped_trimmed_ars_reads.smk"
#include: "rules/sam_to_sortedBam_ars_trimmed.smk"
#include: "rules/featureCounts_quantify_ars_STAR.smk"
#include: "rules/featureCounts_quantify_ars_STAR_softclipped.smk"
#include: "rules/bowtie_map_unmapped_trimmed_ars_reads_multimapping.smk"
#include: "rules/sam_to_sortedBam_ars_trimmed_multimapping.smk"
#include: "rules/featureCounts_quantify_merged_ars_trimmed_multimapping.smk"
##include: "rules/bash_compose_assemble_ars_trimmed_unmapped_reads.smk"
##include: "rules/megahit_assemble_ars_unmapped_reads.smk"
##include: "rules/prodigal_gene_prediction.smk"
##include: "rules/blast_unmapped_reads_against_nr_diamond.smk"
#include: "rules/kraken_taxonomic_identification_unmapped_reads.smk"
