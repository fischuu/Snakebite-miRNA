# Generate more configuration keys
# Define the script_folder dynamically based on the pipeline_folder
SCRIPT_FOLDER = os.path.join(config["pipeline_folder"], "workflow", "scripts")

READS = Path("results/reads/")
WD = os.getcwd()

# references
REF = Path("results/references")
REF_MIRBASE = REF / "mirbase"

# preparation
PREPA = Path("results/preparation")
PREPA_BOWTIE = PREPA / "bowtie"

# preprocess
PRE = Path("results/preprocess/")
CUTADAPT = PRE / "cutadapt/"

# decontamination
DECON = Path("results/decontamination")
DECON_BOWTIE = DECON / "bowtie"
DECON_BOWTIE_TRNA = DECON_BOWTIE / "trna"
DECON_BOWTIE_PHIX = DECON_BOWTIE / "phix"
DECON_SAMTOOLS = DECON / "samtools"
DECON_SAMTOOLS_TRNA = DECON_SAMTOOLS / "trna"
DECON_SAMTOOLS_PHIX = DECON_SAMTOOLS / "phix"

# alignment
ALIGN = Path("results/alignment")
ALIGN_BOWTIE = ALIGN / "bowtie"
ALIGN_BOWTIE_MATURE = ALIGN_BOWTIE / "mature"
ALIGN_BOWTIE_HAIRPIN = ALIGN_BOWTIE / "hairpin"
ALIGN_BOWTIE_MATURE_SPECIES = ALIGN_BOWTIE / "mature_species"
ALIGN_BOWTIE_HAIRPIN_SPECIES = ALIGN_BOWTIE / "hairpin_species"

