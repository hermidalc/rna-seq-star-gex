from os import getcwd
from os.path import join, splitext

DATA_DIR = "data"
RESULTS_DIR = "results"

GENCODE_RESULT_DIR = join(RESULTS_DIR, "gencode")

GENCODE_SPECIES = config["gencode"]["species"]
GENCODE_RELEASE = config["gencode"]["release"]
GENCODE_BUILD = config["gencode"]["build"]

GENCODE_GENOME_FASTA_WRAPPER = join(
    config["wrapper"]["base_url"], "bio/reference/gencode/sequence"
)
GENCODE_GENOME_ANNOT_WRAPPER = join(
    config["wrapper"]["base_url"], "bio/reference/gencode/annotation"
)


rule download_gencode_genome_fasta:
    params:
        species=GENCODE_SPECIES,
        release=GENCODE_RELEASE,
        build=GENCODE_BUILD,
    output:
        GENCODE_GENOME_FASTA_FILE,
    wrapper:
        GENCODE_GENOME_FASTA_WRAPPER


rule download_gencode_genome_annot:
    params:
        species=GENCODE_SPECIES,
        release=GENCODE_RELEASE,
    output:
        GENCODE_GENOME_ANNOT_FILE,
    wrapper:
        GENCODE_GENOME_ANNOT_WRAPPER
