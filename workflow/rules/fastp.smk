
rule run_fastp:
    input:
        fq1=FASTQ1_FILE,
        fq2=FASTQ2_FILE,
    params:
        extra="'--cut_tail --trim_poly_x'",
        # other extra
        # filter reads where 40% of bases have phred quality < 15
        # --unqualified_percent_limit 40
        # filter reads with less than 30% complexity
        # --low_complexity_filter
    output:
        trim_fq1=TRIMMED_FASTQ1_FILE,
        trim_fq2=TRIMMED_FASTQ2_FILE,
        trim_up1=TRIMMED_UNPAIR1_FILE,
        trim_up2=TRIMMED_UNPAIR2_FILE,
        failed=FAILED_READS_FILE,
        html=FASTP_HTML_REPORT_FILE,
        json=FASTP_JSON_REPORT_FILE,
    threads: config["trimming"]["threads"]
    log:
        FASTP_LOG,
    wrapper:
        FASTP_WRAPPER
