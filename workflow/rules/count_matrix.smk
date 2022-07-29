localrules:
    star_count_matrix,


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
