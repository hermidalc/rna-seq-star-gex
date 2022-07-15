def get_fq(wildcards):
    if config["trimming"]["on"]:
        fastq_dir = "data/fastq"
    else:
        fastq_dir = "results/trimmed"
    return {
        "fq1": join(fastq_dir, f"{wildcards.sample}_R1.fastq.gz"),
        "fq2": join(fastq_dir, f"{wildcards.sample}_R2.fastq.gz"),
    }
