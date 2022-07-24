def get_fq(wildcards):
    u = units_df.loc[
        (wildcards.sample, wildcards.unit)
        if hasattr(wildcards, "unit")
        else wildcards.sample,
        ["fq1", "fq2"],
    ].map(
        lambda x: join(
            TRIMMED_RESULTS_DIR if config["trimming"]["activate"] else FASTQ_DIR, x
        )
    )
    return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}


def get_source_fq(wildcards):
    u = units_df.loc[
        (wildcards.sample, wildcards.unit)
        if hasattr(wildcards, "unit")
        else wildcards.sample,
        ["fq1", "fq2"],
    ].map(lambda x: join(FASTQ_DATA_DIR, x))
    return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}
