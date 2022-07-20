from os.path import join

DATA_DIR = "data"
LOG_DIR = "logs"
RESULTS_DIR = "results"
REPORTS_DIR = "reports"

FASTQ_DATA_DIR = join(DATA_DIR, "fastq")
TRIMMED_RESULTS_DIR = join(RESULTS_DIR, config["trimming"]["dirname"])
FASTP_LOG_DIR = join(LOG_DIR, "fastp")
FASTP_REPORTS_DIR = join(REPORTS_DIR, "fastp")

FASTQ1_FILE = join(FASTQ_DIR, "{sample}_R1.fastq.gz")
FASTQ2_FILE = join(FASTQ_DIR, "{sample}_R2.fastq.gz")
TRIMMED_FASTQ1_FILE = join(TRIMMED_RESULTS_DIR, "{sample}_R1.trimmed.fastq.gz")
TRIMMED_FASTQ2_FILE = join(TRIMMED_RESULTS_DIR, "{sample}_R2.trimmed.fastq.gz")
TRIMMED_UNPAIRED_FILE = join(TRIMMED_RESULTS_DIR, "{sample}_UP.trimmed.fastq.gz")
FAILED_READS_FILE = join(TRIMMED_RESULTS_DIR, "{sample}.failed.fastq.gz")
JSON_REPORT_FILE = join(FASTP_REPORTS_DIR, "{sample}_report.json")
HTML_REPORT_FILE = join(FASTP_REPORTS_DIR, "{sample}_report.html")

FASTP_LOG = join(FASTP_LOG_DIR, "{sample}.log")
FASTP_WRAPPER = join(config["wrapper"]["base_url"], "main/bio/fastp")


rule run_fastp:
    input:
        fq1=FASTQ1_FILE,
        fq2=FASTQ2_FILE,
    output:
        trim_fq1=TRIMMED_FASTQ1_FILE,
        trim_fq1=TRIMMED_FASTQ2_FILE,
        trim_up1=TRIMMED_UNPAIRED_FILE,
        trim_up2=TRIMMED_UNPAIRED_FILE,
        failed=FAILED_READS_FILE,
        html=HTML_REPORT_FILE,
        json=JSON_REPORT_FILE,
        extra="'--cut_tail --trim_poly_x'",
        # other extra
        # filter reads where 40% of bases have phred quality < 15
        # --unqualified_percent_limit 40
        # filter reads with less than 30% complexity
        # --low_complexity_filter
    threads: config["trimming"]["threads"]
    log:
        FASTP_LOG,
    wrapper:
        FASTP_WRAPPER
