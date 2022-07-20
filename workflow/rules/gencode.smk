import gzip
from shutil import copyfileobj


rule download_gencode_genome_seq:
    params:
        species=GENCODE_SPECIES,
        release=GENCODE_RELEASE,
        build=GENCODE_BUILD,
        regions=GENCODE_REGIONS,
    output:
        GENCODE_GENOME_SEQ_GZ_FILE,
    wrapper:
        GENCODE_GENOME_SEQ_WRAPPER


rule gunzip_gencode_genome_seq:
    input:
        GENCODE_GENOME_SEQ_GZ_FILE,
    output:
        GENCODE_GENOME_SEQ_FILE,
    run:
        with gzip.open(input, "rb") as f_in:
            with open(output, "wb") as f_out:
                copyfileobj(f_in, f_out)


rule download_gencode_genome_annot:
    params:
        species=GENCODE_SPECIES,
        release=GENCODE_RELEASE,
        regions=GENCODE_REGIONS,
        annot_fmt=GENCODE_ANNOT_FMT,
    output:
        GENCODE_GENOME_ANNOT_GZ_FILE,
    wrapper:
        GENCODE_GENOME_ANNOT_WRAPPER


rule gunzip_gencode_genome_annot:
    input:
        GENCODE_GENOME_ANNOT_GZ_FILE,
    output:
        GENCODE_GENOME_ANNOT_FILE,
    run:
        with gzip.open(input, "rb") as f_in:
            with open(output, "wb") as f_out:
                copyfileobj(f_in, f_out)
