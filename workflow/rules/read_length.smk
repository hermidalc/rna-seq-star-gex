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
        print(f"Getting maximum readlength from {input[0]}")
        df = pd.read_csv(input[0], sep="\t", comment="#", names=cols)
        with open(output[0], "w") as fh:
            fh.write(df["length"].max().astype(str))
        with open(log[0], "w") as fh:
            fh.write("Wrote maximum readlength to file")
