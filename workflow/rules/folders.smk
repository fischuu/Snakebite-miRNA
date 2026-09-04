# Generate more configuration keys
SCRIPT_FOLDER = os.path.join(config["pipeline_folder"], "workflow", "scripts")
WD = os.getcwd()

# Guaranteed-available scratch fallback for tools (e.g. GNU sort -T) that need
# a writable tmp directory: nvme_storage/tmp_storage (config.yaml) may point at
# per-job cluster paths that aren't necessarily visible inside a container, so
# rules that use this probe those first and fall back to this project-local
# directory, which always exists and is on the same storage as results/.
PROJECT_TMP = Path("tmp/")

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
