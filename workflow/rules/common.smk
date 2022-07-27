from os.path import basename


def get_fq(wildcards):
    u = UNIT_DF.loc[
        (wildcards.sample, wildcards.unit)
        if hasattr(wildcards, "unit")
        else wildcards.sample,
        ["fq1", "fq2"],
    ].map(
        lambda x: join(
            TRIMMED_RESULTS_DIR if config["trimming"]["activate"] else FASTQ_PARENT_DIR,
            basename(x),
        )
    )
    return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}


def get_source_fq(wildcards):
    u = UNIT_DF.loc[
        (wildcards.sample, wildcards.unit)
        if hasattr(wildcards, "unit")
        else wildcards.sample,
        ["fq1", "fq2"],
    ].map(lambda x: join(FASTQ_PARENT_DIR, x))
    return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}
