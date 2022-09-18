from os.path import basename, join


def get_fq(wildcards, trimmed):
    u = UNIT_DF.loc[
        (wildcards.sample, wildcards.unit)
        if hasattr(wildcards, "unit")
        else wildcards.sample,
        ["fq1", "fq2"],
    ].map(
        lambda x: join(TRIMMED_RESULTS_DIR, basename(x))
        if config["trim"]["activate"] and trimmed
        else join(FASTQ_DATA_DIR, x)
    )
    return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}
