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
        bam=STAR_BAM_FILE,
    params:
        sample_size=config["rseqc"]["infer_exp"]["sample_size"],
    output:
        infer=RSEQC_INFER_EXPERIMENT_FILE,
        strand=RSEQC_STRAND_INFO_FILE,
    log:
        RSEQC_INFER_EXPERIMENT_LOG,
    wrapper:
        RSEQC_INFER_EXPERIMENT_WRAPPER
