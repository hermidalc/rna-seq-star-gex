import pandas as pd


localrules:
    read_length_histogram,
    max_read_length,


rule read_length_histogram:
    input:
        unpack(get_fq),
    output:
        READ_LENGTH_HISTOGRAM_FILE,
    log:
        READ_LENGTH_HISTOGRAM_LOG,
    wrapper:
        READ_LENGTH_WRAPPER


rule max_read_length:
    input:
        READ_LENGTH_HISTOGRAM_FILE,
    output:
        READ_LENGTH_FILE,
    log:
        READ_LENGTH_LOG,
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
        print(f"Getting maximum read length from {input[0]}", flush=True)
        df = pd.read_csv(input[0], sep="\t", comment="#", names=cols)
        max_length = df["length"].max().astype(str)
        with open(output[0], "w") as fh:
            fh.write(max_length)
        with open(log[0], "w") as fh:
            fh.write(f"Wrote max length {max_length} file")