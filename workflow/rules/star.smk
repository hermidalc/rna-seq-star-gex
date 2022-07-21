import re
from os.path import join


def get_readlength(wildcards):
    with open(expand(READLENGTH_FILE, **wildcards)[0], "r") as fh:
        return re.sub("\D+", "", fh.readline())


localrules:
    run_star_filter_pass1_sj,


rule create_star_genome_index:
    input:
        GENCODE_GENOME_SEQ_FILE,
    output:
        directory(STAR_GENOME_DIR),
    threads: config["star"]["index"]["threads"]
    log:
        STAR_GENOME_LOG,
    wrapper:
        STAR_GENOME_WRAPPER


rule run_star_align_pass1:
    input:
        unpack(get_fq),
        index=STAR_GENOME_DIR,
        gtf=GENCODE_GENOME_ANNOT_FILE,
    params:
        out_dir=STAR_PASS1_OUTPUT_DIR,
        readlength=get_readlength,
        extra=f"--outSAMtype None",
    output:
        STAR_PASS1_SJ_FILE,
    threads: config["star"]["index"]["threads"]
    log:
        STAR_ALIGN_PASS1_LOG,
    wrapper:
        STAR_ALIGN_WRAPPER


rule run_star_filter_pass1_sj:
    input:
        STAR_PASS1_SJ_FILE,
    output:
        STAR_PASS1_SJ_FILTERED_FILE,
    log:
        STAR_PASS1_SJ_FILTERED_LOG,
    run:
        num_novel_sj = 0
        num_filtered_novel_sj = 0
        # chromosomal and non-mitochondrial (regex specific to GTF style!)
        chr_no_mt_regex = re.compile("chr([1-9][0-9]?|X|Y)")
        with open(input[0], "r") as f_in:
            with open(output[0], "w") as f_out:
                for line in f_in:
                    fields = line.rstrip().split("\s+")
                    # skip annotated (and keep novel, since already get added from GTF)
                    if fields[5] != 0:
                        continue
                    num_novel_sj += 1
                    if (
                        # chromosomal and non-mitochondrial
                        re.match(chr_no_mt_regex, fields[0])
                        and
                        # canonical
                        fields[4] > 0
                        and
                        # supported by at least one unique mapper
                        fields[6] > 0
                    ):
                        f_out.write(line)
                        num_filtered_novel_sj += 1
        with open(log[0], "w") as fh:
            fh.write(f"Num novel sj: {num_novel_sj:d}\n")
            fh.write(f"Num filtered novel sj: {num_filtered_novel_sj:d}\n")


rule run_star_align_pass2:
    input:
        unpack(get_fq),
        index=STAR_GENOME_DIR,
        gtf=GENCODE_GENOME_ANNOT_FILE,
        sjdb=STAR_PASS1_SJ_FILTERED_FILE,
    params:
        out_dir=STAR_PASS2_OUTPUT_DIR,
        readlength=get_readlength,
        extra=(
            f"--outSAMattrRGline {SAM_ATTR_RG_LINE}"
            " --outFilterType BySJout"
            " --outSAMattributes All"
            " --outSAMstrandField intronMotif"
            f" --outSAMtype BAM {STAR_BAM_SORT}"
            " --outSAMunmapped Within"
            " --quantMode ReadCounts"
        ),
    output:
        bam_file=STAR_PASS2_BAM_FILE,
        count_file=STAR_PASS2_READCOUNT_FILE,
    log:
        STAR_ALIGN_PASS2_LOG,
    threads: config["star"]["align"]["threads"]
    wrapper:
        STAR_ALIGN_WRAPPER
