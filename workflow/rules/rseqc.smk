rule gffutils_gtf2bed:
    input:
        GENCODE_GENOME_ANNOT_FILE,
    output:
        bed=RSEQC_GENOME_ANNOT_FILE,
        db=temp(RSEQC_GENOME_ANNOT_DB),
    log:
        RSEQC_ANNOT_LOG,
    wrapper:
        RSEQC_ANNOT_WRAPPER


rule rseqc_infer_experiment:
    input:
        bed=RSEQC_GENOME_ANNOT_FILE,
        bam=STAR_PASS2_BAM_FILE,
    output:
        RSEQC_INFER_EXPERIMENT_FILE,
    log:
        RSEQC_INFER_EXPERIMENT_LOG,
    priority: 1
    wrapper:
        RSEQC_INFER_EXPERIMENT_WRAPPER
