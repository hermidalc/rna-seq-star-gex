rule gsea_msigdb:
    conda:
        "../envs/gsea_msigdb.yaml"
    input:
        DIFFEXP_RESULTS_TABLE_FILE,
    params:
        de_meth="{de_meth}",
        species=GENCODE_SPECIES,
        seed=config["random_seed"],
        gsea_padj=config["gsea"]["padj"],
        db_cat=lambda wc: next(
            db["category"]
            for db in config["gsea"]["msigdb"]
            if wc.msigdb
            == (
                db["category"].lower()
                if db["subcategory"] is None
                else (
                    db["category"].lower()
                    + "_"
                    + db["subcategory"].lower().replace(":", "_")
                )
            )
        ),
        db_sub=lambda wc: next(
            db["subcategory"]
            for db in config["gsea"]["msigdb"]
            if wc.msigdb
            == (
                db["category"].lower()
                if db["subcategory"] is None
                else (
                    db["category"].lower()
                    + "_"
                    + db["subcategory"].lower().replace(":", "_")
                )
            )
        ),
        collapse=lambda wc: next(
            db["collapse"]
            for db in config["gsea"]["msigdb"]
            if wc.msigdb
            == (
                db["category"].lower()
                if db["subcategory"] is None
                else (
                    db["category"].lower()
                    + "_"
                    + db["subcategory"].lower().replace(":", "_")
                )
            )
        ),
    output:
        GSEA_MSIGDB_RESULTS_FILE,
    log:
        READ_LENGTH_LOG,
    script:
        "../scripts/gsea_msigdb.R"
