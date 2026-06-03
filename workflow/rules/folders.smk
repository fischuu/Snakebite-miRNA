# Generate more configuration keys
# Define the script_folder dynamically based on the pipeline_folder
SCRIPT_FOLDER = os.path.join(config["pipeline_folder"], "workflow", "scripts")

READS = Path("results/reads/")
WD = os.getcwd()

# reference
REFERENCE = Path("results/reference/")
#HOSTS = REFERENCE / "hosts"

# preprocess
PRE = Path("results/preprocess/")
CUTADAPT = PRE / "cutadapt/"
