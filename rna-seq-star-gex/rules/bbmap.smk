from os.path import join

import pandas as pd

BBMAP_DATA_DIR = join(DATA_DIR, "bbmap")
FASTQ_DIR = join(DATA_DIR, "fastq")

READ_LENGTH_HISTOGRAM = join(BBMAP_DATA_DIR, "read_length_histogram.tsv")


localrules:
    get_readlength_histogram,
    get_max_readlength,


rule get_readlength_histogram:
    input:
        unpack(get_fq),
    output:
        READ_LENGTH_HISTOGRAM,
    wrapper:
        "https://github.com/hermidalc/snakemake-wrappers/tree/main/bio/bbmap/readlength"


rule get_max_readlength:
    conda:
        "../envs/bbmap.yaml"
    input:
        READ_LENGTH_HISTOGRAM,
    output:
        READ_LENGTH_FILE,
    run:
        cols = [
            "Length",
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
        max_read_length = df["Length"].max().astype(int)
        with open(output, "wb") as fh:
            print(fh, max_read_length, end="")
