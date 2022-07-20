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
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


def get_fq(wildcards):
    if config["trimming"]["activate"]:
        fastq_dir = TRIMMED_RESULTS_DIR
    else:
        fastq_dir = FASTQ_DATA_DIR
    return {
        zip(
            ["fq1", "fq2"],
            expand(
                join(fastq_dir, f"{FASTQ_PREFIX_WILDCARD_STR}_{{pair}}.fastq.gz"),
                pair=FASTQ_PAIRS,
                **wildcards,
            ),
        )
    }
