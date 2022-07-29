localrules:
    count_eset,


rule count_eset:
    input:
        assay=COUNT_MATRIX_FILE,
        pheno=SAMPLE_CONFIG_FILE,
        annot=GENCODE_GENE_ANNOT_FILE,
    params:
        samples=SAMPLE_LABELS,
    output:
        COUNT_ESET_FILE,
    log:
        COUNT_ESET_LOG,
    wrapper:
        COUNT_ESET_WRAPPER
