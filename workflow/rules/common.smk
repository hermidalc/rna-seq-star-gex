import re

import pandas as pd
from snakemake.utils import validate

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)
validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(
        ["sample_name", "unit_name"]
        if re.search("\{sample\}", FASTQ_PREFIX_WILDCARD_STR)
        else "unit_name",
        drop=False,
    )
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


def get_fq(wildcards):
    u = units.loc[
        (wildcards.sample, wildcards.unit)
        if hasattr(wildcards, "sample")
        else wildcards.unit,
        ["fq1", "fq2"],
    ].map(
        lambda x: join(
            TRIMMED_RESULTS_DIR if config["trimming"]["activate"] else FASTQ_DIR, x
        )
    )
    return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}


def get_source_fq(wildcards):
    u = units.loc[
        (wildcards.sample, wildcards.unit)
        if hasattr(wildcards, "sample")
        else wildcards.unit,
        ["fq1", "fq2"],
    ].map(lambda x: join(FASTQ_DATA_DIR, x))
    return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}
