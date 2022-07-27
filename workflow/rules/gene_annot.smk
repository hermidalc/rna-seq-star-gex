rule gtf_gene_annot:
    input:
        GENCODE_GENOME_ANNOT_FILE,
    params:
        length_col=config["gencode"]["gene_annot"]["length_col"],
    output:
        GENCODE_GENE_ANNOT_FILE,
    log:
        GENCODE_GENE_ANNOT_LOG,
    wrapper:
        GENCODE_GENE_ANNOT_WRAPPER
