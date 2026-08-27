# Generate more configuration keys
SCRIPT_FOLDER = os.path.join(config["pipeline_folder"], "workflow", "scripts")
WD = os.getcwd()

# reads
READS = Path("results/reads/")

# reference
REFERENCE = Path("results/reference/")
MIRBASE = REFERENCE / "mirbase/"
STAR_INDEX = REFERENCE / "star_index/"
STAR_INDEX_MATURE = STAR_INDEX / "mature"
STAR_INDEX_MATURE_SPECIES = STAR_INDEX / "mature_species"
STAR_INDEX_HAIRPIN = STAR_INDEX / "hairpin"
STAR_INDEX_HAIRPIN_SPECIES = STAR_INDEX / "hairpin_species"
STAR_INDEX_GENOME = STAR_INDEX / "genome"
STAR_INDEX_TMP = STAR_INDEX / "tmp/"

# preprocess
PRE = Path("results/preprocess/")
TRIMMED = PRE / "cutadapt/"
CONCATENATED = PRE / "concatenate/"

# decontaminate
DECONTAM = Path("results/decontaminate/")
TRNA = DECONTAM / "trna/"
PHIX = DECONTAM / "phix/"

# align
ALIGN = Path("results/align/")
MATURE_BOWTIE = ALIGN / "bowtie/mature/"
MATURE_SPECIES_BOWTIE = ALIGN / "bowtie/mature_species/"
HAIRPIN_BOWTIE = ALIGN / "bowtie/hairpin/"
GENOME_BOWTIE = ALIGN / "bowtie/genome/"
MATURE_STAR = ALIGN / "star/mature/"
MATURE_SPECIES_STAR = ALIGN / "star/mature_species/"
HAIRPIN_STAR = ALIGN / "star/hairpin/"
HAIRPIN_SPECIES_STAR = ALIGN / "star/hairpin_species/"
GENOME_STAR = ALIGN / "star/genome/"
STAR_ALIGN_TMP = ALIGN / "star_tmp/"
SOFTCLIPPED = ALIGN / "softclipped/"

# novel_mirna
NOVEL_MIRNA = Path("results/novel_mirna/")
MPILEUP = NOVEL_MIRNA / "mpileup/"
NOVEL_LOCI = NOVEL_MIRNA / "bedtools/"

# quantify
QUANT = Path("results/quantify/")
QUANT_BOWTIE = QUANT / "bowtie/"
QUANT_STAR = QUANT / "star/"

# reporting
PIPELINE_REPORT = Path("Rreports/")
