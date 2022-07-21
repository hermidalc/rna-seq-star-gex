rule download_gencode_genome_seq:
    params:
        species=GENCODE_SPECIES,
        release=GENCODE_RELEASE,
        build=GENCODE_BUILD,
        regions=GENCODE_REGIONS,
    output:
        GENCODE_GENOME_SEQ_FILE,
    wrapper:
        GENCODE_GENOME_SEQ_WRAPPER


rule download_gencode_genome_annot:
    params:
        species=GENCODE_SPECIES,
        release=GENCODE_RELEASE,
        regions=GENCODE_REGIONS,
        annot_fmt=GENCODE_ANNOT_FMT,
    output:
        GENCODE_GENOME_ANNOT_FILE,
    wrapper:
        GENCODE_GENOME_ANNOT_WRAPPER
