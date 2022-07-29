from os.path import basename


def get_fq(wildcards, trimmed):
    u = UNIT_DF.loc[
        (wildcards.sample, wildcards.unit)
        if hasattr(wildcards, "unit")
        else wildcards.sample,
        ["fq1", "fq2"],
    ].map(
        lambda x: join(
            TRIMMED_RESULTS_DIR
            if config["trimming"]["activate"] and trimmed
            else FASTQ_PARENT_DIR,
            basename(x),
        )
    )
    return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}
