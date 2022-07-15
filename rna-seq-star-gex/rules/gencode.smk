from os import getcwd
from os.path import join, splitext

WORK_DIR = getcwd()
DATA_DIR = join(WORK_DIR, "data")
GENCODE_DIR = join(DATA_DIR, "gencode")

GENCODE_BASE_URL = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/"
GENCODE_MOUSE_VERSION = "M30"
GENCODE_MOUSE_BUILD = "GRCm39"

GENCODE_MOUSE_BASE_URL = join(
    GENCODE_BASE_URL, f"Gencode_mouse/release_{MOUSE_GENCODE_VERSION}"
)
GENCODE_MOUSE_GENOME_FASTA_GZ_URL = join(
    GENCODE_MOUSE_BASE_URL, GENCODE_MOUSE_GENOME_FASTA_GZ_FILENAME
)

GENCODE_MOUSE_GTF_GZ_URL = join(
    GENCODE_MOUSE_BASE_URL, GENCODE_MOUSE_GENOME_GTF_GZ_FILENAME
)

GENCODE_MOUSE_GENOME_FASTA_FILE = join(
    GENCODE_DIR, splitext(GENCODE_MOUSE_GENOME_FASTA_GZ_FILENAME)[0]
)
GENCODE_MOUSE_GTF_FILE = join(
    GENCODE_DIR, splitext(GENCODE_MOUSE_GENOME_GTF_GZ_FILENAME)[0]
)


rule download_gencode_genome_fasta:
    params:
        GENCODE_MOUSE_GENOME_FASTA_GZ_URL,
    output:
        GENCODE_MOUSE_GENOME_FASTA_FILE,
    shell:
        "wget -O - {params} | gunzip -c > {output}"


rule download_gencode_gtf:
    params:
        GENCODE_MOUSE_GTF_GZ_URL,
    output:
        GENCODE_MOUSE_GTF_FILE,
    shell:
        "wget -O - {params} | gunzip -c > {output}"
