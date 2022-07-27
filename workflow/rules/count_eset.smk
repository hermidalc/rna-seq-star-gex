localrules:
    count_eset,


rule count_eset:
    input:
        assay=COUNT_MATRIX_FILE,
        pheno=SAMPLE_DF,
        feature=GENCODE_GENE_ANNOT_FILE,
    output:
        COUNT_ESET_FILE,
    log:
        COUNT_ESET_LOG,
    wrapper:
        COUNT_ESET_WRAPPER
