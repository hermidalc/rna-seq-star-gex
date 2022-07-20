import pandas as pd
import snakemake

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
with open(snakemake.output, "wb") as fh:
    print(fh, max_readlength, end="")
