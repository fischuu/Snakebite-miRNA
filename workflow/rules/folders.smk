# Generate more configuration keys
SCRIPT_FOLDER = os.path.join(config["pipeline_folder"], "workflow", "scripts")
WD = os.getcwd()

RAW = Path("FASTQ/RAW/")

# reference
REFDIR = Path("References/")
STARBASE = Path("References/STAR/")
STAR_INDEX_MATURE = STARBASE / "Mature"
STAR_INDEX_MATURE_SPECIES = STARBASE / "MatureSpecies"
STAR_INDEX_HAIRPIN = STARBASE / "Hairpin"
STAR_INDEX_HAIRPIN_SPECIES = STARBASE / "HairpinSpecies"
STAR_INDEX_GENOME = STARBASE / "Reference"
STAR_TMP = Path("STAR_tmp/")

# preprocess
TRIMMED = Path("FASTQ/TRIMMED/")
CONCATENATED = Path("FASTQ/CONCATENATED/")

# decontaminate
TRNA_FASTQ = Path("FASTQ/tRNA/")
TRNA_BAM = Path("BAM/BOWTIE/tRNA/")
TRNA_STATS = Path("STATS/BOWTIE/tRNA/")
PHIX_FASTQ = Path("FASTQ/PhiX/")
PHIX_BAM = Path("BAM/BOWTIE/PhiX/")
PHIX_STATS = Path("STATS/BOWTIE/PhiX/")

# align
MATURE_FASTQ = Path("FASTQ/Mature/")
MATURE_BAM_BOWTIE = Path("BAM/BOWTIE/Mature/")
MATURE_BAM_STAR = Path("BAM/STAR/Mature/")
MATURE_STATS_BOWTIE = Path("STATS/BOWTIE/Mature/")
MATURE_STATS_STAR = Path("STATS/STAR/Mature/")

MATURE_SPECIES_FASTQ = Path("FASTQ/Mature_Species/")
MATURE_SPECIES_BAM_BOWTIE = Path("BAM/BOWTIE/Mature_Species/")
MATURE_SPECIES_BAM_STAR = Path("BAM/STAR/Mature_Species/")
MATURE_SPECIES_STATS_BOWTIE = Path("STATS/BOWTIE/Mature_Species/")
MATURE_SPECIES_STATS_STAR = Path("STATS/STAR/Mature_Species/")

HAIRPIN_FASTQ = Path("FASTQ/Hairpin/")
HAIRPIN_BAM_BOWTIE = Path("BAM/BOWTIE/Hairpin/")
HAIRPIN_BAM_STAR = Path("BAM/STAR/Hairpin/")
HAIRPIN_STATS_BOWTIE = Path("STATS/BOWTIE/Hairpin/")
HAIRPIN_STATS_STAR = Path("STATS/STAR/Hairpin/")

HAIRPIN_SPECIES_FASTQ = Path("FASTQ/Hairpin_Species/")
HAIRPIN_SPECIES_BAM_STAR = Path("BAM/STAR/Hairpin_Species/")
HAIRPIN_SPECIES_STATS_STAR = Path("STATS/STAR/Hairpin_Species/")

GENOME_FASTQ = Path("FASTQ/Reference/")
GENOME_BAM_BOWTIE = Path("BAM/BOWTIE/Reference/")
GENOME_BAM_STAR = Path("BAM/STAR/Reference/")
GENOME_STATS_BOWTIE = Path("STATS/BOWTIE/Reference/")
GENOME_STATS_STAR = Path("STATS/STAR/Reference/")

SOFTCLIPPED = Path("FASTA/STAR/")

# quantify
QUANT_BOWTIE = Path("QUANTIFICATION/BOWTIE/")
QUANT_STAR = Path("QUANTIFICATION/STAR/")

# novel_mirna
MPILEUP = Path("MPILEUP/mpileup_reference/")

# QC
QC = Path("QC/")

# reporting
PIPELINE_REPORT = Path("PipelineReport/")
