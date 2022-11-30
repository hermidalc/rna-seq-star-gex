rule star_genome_index:
    input:
        fastas=GENCODE_GENOME_SEQ_FILE,
    output:
        directory(STAR_GENOME_DIR),
    resources:
        tmpdir=TEMP_DIR,
    threads: STAR_INDEX_THREADS
    log:
        STAR_GENOME_LOG,
    wrapper:
        STAR_GENOME_WRAPPER


rule star_align_pass1:
    input:
        unpack(lambda w: get_fq(w, trimmed=True)),
        gtf=GENCODE_GENOME_ANNOT_FILE,
        index=STAR_GENOME_DIR,
        read_length=READ_LENGTH_FILE,
    params:
        out_dir=STAR_PASS1_OUTPUT_DIR,
        extra=(
            " --alignIntronMax 1000000"
            " --alignIntronMin 20"
            " --alignMatesGapMax 1000000"
            " --alignSJDBoverhangMin 1"
            " --alignSJoverhangMin 8"
            " --alignSoftClipAtReferenceEnds Yes"
            " --limitSjdbInsertNsj 1200000"
            " --outFilterIntronMotifs None"
            " --outFilterMismatchNmax 999"
            " --outFilterMismatchNoverLmax 0.1"
            " --outFilterMultimapNmax 20"
            " --outSAMtype None"
        ),
    output:
        STAR_PASS1_SJ_FILE,
    resources:
        tmpdir=TEMP_DIR,
    threads: STAR_ALIGN_THREADS
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
        unpack(lambda w: get_fq(w, trimmed=True)),
        index=STAR_GENOME_DIR,
        gtf=GENCODE_GENOME_ANNOT_FILE,
        read_length=READ_LENGTH_FILE,
        sj=STAR_PASS1_SJ_FILTERED_FILE,
    params:
        out_dir=STAR_PASS2_OUTPUT_DIR,
        extra=(
            " --alignIntronMax 1000000"
            " --alignIntronMin 20"
            " --alignMatesGapMax 1000000"
            " --alignSJDBoverhangMin 1"
            " --alignSJoverhangMin 8"
            " --alignSoftClipAtReferenceEnds Yes"
            " --chimJunctionOverhangMin 15"
            " --chimMainSegmentMultNmax 1"
            " --chimSegmentMin 15"
            " --chimOutType Junctions SeparateSAMold WithinBAM SoftClip"
            " --limitSjdbInsertNsj 1200000"
            " --outFilterIntronMotifs None"
            " --outFilterMismatchNmax 999"
            " --outFilterMismatchNoverLmax 0.1"
            " --outFilterMultimapNmax 20"
            " --outFilterType BySJout"
            f" --outSAMattrRGline {STAR_SAM_ATTR_RG_LINE}"
            f" --outSAMtype BAM {STAR_BAM_SORT}"
            " --outSAMstrandField intronMotif"
            " --outSAMattributes All"
            " --outSAMunmapped Within"
            " --quantMode GeneCounts"
        ),
    output:
        bam_file=STAR_BAM_FILE,
        count_file=STAR_READ_COUNT_FILE,
    log:
        STAR_ALIGN_PASS2_LOG,
    resources:
        tmpdir=TEMP_DIR,
    threads: STAR_ALIGN_THREADS
    wrapper:
        STAR_ALIGN_WRAPPER
