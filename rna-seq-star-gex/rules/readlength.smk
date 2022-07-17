from os.path import join

import pandas as pd

LOG_DIR = "logs"
RESULTS_DIR = "results"

READLENGTH_LOG_DIR = join(LOG_DIR, "readlength")
READLENGTH_RESULTS_DIR = join(RESULTS_DIR, "readlength")

READLENGTH_HISTOGRAM_FILE = join(READLENGTH_RESULTS_DIR, "{sample}_histogram.tsv")
READLENGTH_FILE = join(READLENGTH_RESULTS_DIR, "{sample}_length.txt")

READLENGTH_HISTOGRAM_LOG = join(READLENGTH_LOG_DIR, "{sample}_histogram.log")
READLENGTH_LOG = join(READLENGTH_LOG_DIR, "{sample}_length.log")


localrules:
    get_readlength_histogram,
    get_max_readlength,


rule get_readlength_histogram:
    input:
        unpack(get_fq),
    output:
        READLENGTH_HISTOGRAM_FILE,
    log:
        READLENGTH_HISTOGRAM_LOG,
    wrapper:
        "https://github.com/hermidalc/snakemake-wrappers/tree/main/bio/bbmap/readlength"


rule get_max_readlength:
    conda:
        "../envs/readlength.yaml"
    input:
        READLENGTH_HISTOGRAM_FILE,
    output:
        READLENGTH_FILE,
    log:
        READLENGTH_LOG,
    run:
        cols = [
            "length",
            "reads",
            "pct_reads",
            "cum_reads",
            "cum_pct_reads",
            "bases",
            "pct_bases",
            "cum_bases",
            "cum_pct_bases",
        ]
        df = pd.read_csv(input, sep="\t", comment="#", names=cols)
        max_readlength = df["length"].max().astype(int)
        with open(output, "wb") as fh:
            print(fh, max_readlength, end="")
