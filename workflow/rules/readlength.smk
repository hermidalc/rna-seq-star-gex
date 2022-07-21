import pandas as pd


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
        READLENGTH_WRAPPER


rule get_max_readlength:
    input:
        READLENGTH_HISTOGRAM_FILE,
    output:
        READLENGTH_FILE,
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
        df = pd.read_csv(input[0], sep="\t", comment="#", names=cols)
        with open(output[0], "w") as fh:
            fh.write(df["length"].max().astype(int))
