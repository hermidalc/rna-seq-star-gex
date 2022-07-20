import pandas as pd


localrules:
    get_readlength_histogram,
    get_max_readlength,


rule get_readlength_histogram:
    input:
        unpack(get_fq),
    output:
        READLENGTH_HISTOGRAM_FILE,
    log:
        READLENGTH_HISTOGRAM_LOG,
    wrapper:
        READLENGTH_WRAPPER


rule get_max_readlength:
    conda:
        "../envs/readlength.yaml"
    input:
        READLENGTH_HISTOGRAM_FILE,
    output:
        READLENGTH_FILE,
    log:
        READLENGTH_LOG,
    script:
        "../scripts/get_max_readlength.py"
