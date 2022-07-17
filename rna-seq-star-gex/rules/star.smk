import re
from os.path import join

DATA_DIR = "data"
LOG_DIR = "logs"
RESULTS_DIR = "results"

GTF_FILE = join(DATA_DIR, config["ref"]["gtf"])

STAR_LOG_DIR = join(LOG_DIR, "star")
STAR_RESULTS_DIR = join(RESULTS_DIR, "star")
STAR_GENOME_DIR = join(STAR_RESULTS_DIR, "{genome}")
STAR_OUTPUT_DIR = join(STAR_RESULTS_DIR, "{sample}")
STAR_PASS1_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARpass1")
STAR_PASS2_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARpass2")

STAR_PASS1_SJ_FILE = join(STAR_PASS1_OUTPUT_DIR, "SJ.out.tab")
STAR_PASS1_SJ_FILTERED_FILE = join(STAR_PASS1_OUTPUT_DIR, "SJ.filtered.out.tab")

STAR_PASS2_BAM_FILE = join(STAR_PASS2_OUTPUT_DIR, "Aligned.sortedByCoord.out.bam")
STAR_PASS2_READCOUNT_FILE = join(STAR_PASS2_OUTPUT_DIR, "ReadsPerGene.out.tab")
STAR_BAM_FILE = join(STAR_PE_OUTPUT_DIR, "Aligned.sortedByCoord.out.bam")
STAR_READCOUNT_FILE = join(STAR_PE_OUTPUT_DIR, "ReadsPerGene.out.tab")

SAM_ATTR_RG_LINE = "--outSAMattrRGline ID:{sample} PL:Illumina SM:{sample} LB:RNA"

STAR_GENOME_LOG = join(STAR_LOG_DIR, "{genome}.log")
STAR_ALIGN_LOG = join(STAR_LOG_DIR, "{sample}.log")


def get_readlength(wildcards):
    file = join(READLENGTH_RESULTS_DIR, f"{wildcards.sample}_length.txt")
    with open(file, "rb") as fh:
        return fh.readline().strip().replace(" ", "")


localrules:
    filter_star_sj_pass1,


rule create_star_genome_index:
    input:
        STAR_GENOME_FASTA,
    output:
        directory(STAR_GENOME_DIR),
    threads: config["star"]["index"]["threads"]
    log:
        STAR_GENOME_LOG,
    wrapper:
        "https://github.com/hermidalc/snakemake-wrappers/tree/main/bio/star/index"


rule run_star_pe_pass1:
    input:
        unpack(get_fq),
        index=STAR_GENOME_DIR,
        gtf=GTF_FILE,
        extra=f"{SAM_ATTR_RG_LINE} --outSAMtype None",
    output:
        STAR_PASS1_SJ_FILE,
    params:
        readlength=get_readlength,
    threads: config["star"]["index"]["threads"]
    log:
        STAR_ALIGN_LOG,
    wrapper:
        "https://github.com/hermidalc/snakemake-wrappers/tree/main/bio/star/align"


rule run_star_filter_pass1_sj:
    input:
        STAR_PASS1_SJ_FILE,
    output:
        STAR_PASS1_SJ_FILTERED_FILE,
    run:
        num_novel_sj = 0
        num_filtered_novel_sj = 0
        # chromosomal and non-mitochondrial (regex specific to GTF style!)
        chr_no_mt_regex = re.compile("chr([1-9][0-9]?|X|Y)")
        with open(input, "rb") as f_in:
            with open(output, "wb") as f_out:
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
                        print(f_out, line)
                        num_filtered_novel_sj += 1


rule run_star_pe_pass2:
    input:
        unpack(get_fq),
        index=STAR_GENOME_DIR,
        gtf=GTF_FILE,
        sjdb=STAR_PASS1_SJ_FILTERED_FILE,
        extra=(
            f"--outSAMattrRGline {SAM_ATTR_RG_LINE} "
            "--outFilterType BySJout "
            "--outSAMattributes NH HI AS nM NM ch "
            "--outSAMstrandField intronMotif "
            "--outSAMtype BAM Unsorted "
            "--outSAMunmapped Within "
            "--quantMode ReadCounts "
        ),
    output:
        bam_file=STAR_PASS2_BAM_FILE,
        count_file=STAR_PASS2_READCOUNT_FILE,
    params:
        readlength=get_readlength,
    threads: config["star"]["align"]["threads"]
    wrapper:
        "https://github.com/hermidalc/snakemake-wrappers/tree/main/bio/star/align"
