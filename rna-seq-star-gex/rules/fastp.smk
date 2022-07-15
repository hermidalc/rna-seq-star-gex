import re
from os import getcwd
from os.path import basename, join

WORK_DIR = getcwd()
DATA_DIR = join(WORK_DIR, "data")
FASTP_DATA_DIR = join(DATA_DIR, "fastp")
FASTQ_DATA_DIR = join(DATA_DIR, "fastq")

FASTQ1_FILE = join(FASTQ_DIR, "{sample}_R1.fastq.gz")
FASTQ2_FILE = join(FASTQ_DIR, "{sample}_R2.fastq.gz")
TRIMMED_FASTQ1_FILE = join(FASTQ_DATA_DIR, "{sample}_R1.trimmed.fastq.gz")
TRIMMED_FASTQ2_FILE = join(FASTQ_DATA_DIR, "{sample}_R2.trimmed.fastq.gz")
TRIMMED_UNPAIRED_FILE = join(FASTQ_DATA_DIR, "{sample}.trimmed.{fastq_ext}")
FAILED_READS_FILE = join(FASTP_DATA_DIR, "{sample}.failed.{fastq_ext}")
FASTP_JSON_REPORT = join(FASTP_DATA_DIR, "{sample}_report.json")
FASTP_HTML_REPORT = join(FASTP_DATA_DIR, "{sample}_report.html")


rule run_fastp:
    input:
        fq1_file=FASTQ1_FILE,
        fq2_file=FASTQ2_FILE,
    output:
        trim_fq1_file=TRIMMED_FASTQ1_FILE,
        trim_fq1_file=TRIMMED_FASTQ2_FILE,
        trim_up_file=TRIMMED_UNPAIRED_FILE,
        failed_file=FAILED_READS_FILE,
        json_report_file=JSON_REPORT_FILE,
        html_report_file=HTML_REPORT_FILE,
    threads: 8
    shell:
        """
        fastp \
        --threads {threads} \
        --in1 {input.fq1_file} \
        --in2 {input.fq2_file} \
        --out1 {output.trim_fq1_file} \
        --out2 {output.trim_fq2_file} \
        --unpaired1 {output.trim_up_file} \
        --unpaired2 {output.trim_up_file} \
        --failed_out {output.failed_file} \
        --json {output.json_report_file} \
        --html {output.html_report_file} \
        --cut_tail \
        --trim_poly_x
        """
        # unused
        # filter reads where 40% of bases have phred quality < 15
        # --unqualified_percent_limit 40
        # filter reads with less than 30% complexity
        # --low_complexity_filter
