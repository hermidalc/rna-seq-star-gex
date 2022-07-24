localrules:
    star_filter_pass1_sj,


rule star_genome_index:
    input:
        fastas=GENCODE_GENOME_SEQ_FILE,
    output:
        directory(STAR_GENOME_DIR),
    threads: config["star"]["index"]["threads"]
    log:
        STAR_GENOME_LOG,
    wrapper:
        STAR_GENOME_WRAPPER


rule star_align_pass1:
    input:
        unpack(get_fq),
        index=STAR_GENOME_DIR,
        read_length=READ_LENGTH_FILE,
    params:
        out_dir=STAR_PASS1_OUTPUT_DIR,
        extra=f"--outSAMtype None",
    output:
        STAR_PASS1_SJ_FILE,
    threads: config["star"]["index"]["threads"]
    log:
        STAR_ALIGN_PASS1_LOG,
    wrapper:
        STAR_ALIGN_WRAPPER


rule star_filter_pass1_sj:
    input:
        STAR_PASS1_SJ_FILE,
    output:
        STAR_PASS1_SJ_FILTERED_FILE,
    log:
        STAR_PASS1_SJ_FILTERED_LOG,
    run:
        print(f"Filtering STAR {input[0]}", flush=True)
        num_novel_sj = 0
        num_filtered_novel_sj = 0
        ## chromosomal and/or non-mitochondrial
        chr_regex = re.compile(config["star"]["sj_filter_chr_regex"])
        with open(input[0], "r") as f_in:
            with open(output[0], "w") as f_out:
                for line in f_in:
                    fields = line.rstrip().split("\t")
                    num_fields = len(fields)
                    assert (
                        num_fields == 9
                    ), "Number of fields {num_fields} not equal to 9"
                    # skip annotated (and keep novel, since already get added from GTF)
                    if int(fields[5]) != 0:
                        continue
                    num_novel_sj += 1
                    if (
                        ## chromosomal and/or non-mitochondrial
                        re.match(chr_regex, fields[0])
                        and
                        # canonical
                        int(fields[4]) > 0
                        and
                        # supported by at least one unique mapper
                        int(fields[6]) > 0
                    ):
                        f_out.write(line)
                        num_filtered_novel_sj += 1
        with open(log[0], "w") as fh:
            fh.write(f"Num novel sj: {num_novel_sj}\n")
            fh.write(f"Num filtered novel sj: {num_filtered_novel_sj}\n")


rule star_align_pass2:
    input:
        unpack(get_fq),
        index=STAR_GENOME_DIR,
        gtf=GENCODE_GENOME_ANNOT_FILE,
        read_length=READ_LENGTH_FILE,
        sj=STAR_PASS1_SJ_FILTERED_FILE,
    params:
        out_dir=STAR_PASS2_OUTPUT_DIR,
        extra=(
            "--chimOutType Junctions SeparateSAMold WithinBAM SoftClip"
            f" --outSAMattrRGline {SAM_ATTR_RG_LINE}"
            f" --outSAMtype BAM {STAR_BAM_SORT}"
            " --outSAMattributes All"
            " --outSAMunmapped Within"
            " --outFilterType BySJout"
            " --quantMode GeneCounts"
        ),
    output:
        bam_file=STAR_BAM_FILE,
        count_file=STAR_READ_COUNT_FILE,
    log:
        STAR_ALIGN_PASS2_LOG,
    threads: config["star"]["align"]["threads"]
    wrapper:
        STAR_ALIGN_WRAPPER
