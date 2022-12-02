rule star_count_matrix:
    input:
        counts=expand(STAR_READ_COUNT_FILE, zip, **EXPAND_PARAMS),
        strand=expand(RSEQC_STRAND_INFO_FILE, zip, **EXPAND_PARAMS),
    params:
        samples=SAMPLE_LABELS,
    output:
        COUNT_MATRIX_FILE,
    log:
        COUNT_MATRIX_LOG,
    wrapper:
        COUNT_MATRIX_WRAPPER


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
