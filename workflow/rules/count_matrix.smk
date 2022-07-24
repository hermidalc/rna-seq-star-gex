localrules:
    star_count_matrix,


rule star_count_matrix:
    input:
        counts=STAR_READ_COUNT_FILE,
        strand=RSEQC_STRAND_INFO_FILE,
    params:
        samples=SAMPLES,
    output:
        COUNT_MATRIX_FILE,
    log:
        COUNT_MATRIX_LOG,
    wrapper:
        COUNT_MATRIX_WRAPPER
