def get_fq(wildcards):
    if config["trimming"]["activate"]:
        fastq_dir = TRIMMED_RESULTS_DIR
    else:
        fastq_dir = FASTQ_DATA_DIR
    return {
        "fq1": join(fastq_dir, f"{wildcards.sample}_R1.fastq.gz"),
        "fq2": join(fastq_dir, f"{wildcards.sample}_R2.fastq.gz"),
    }
