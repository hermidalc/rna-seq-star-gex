from os.path import join

DATA_DIR = "data"
RESULTS_DIR = "results"

STAR_RESULTS_DIR = join(DATA_DIR, "star")
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


localrules:
    filter_star_sj_pass1,


rule create_star_genome_index:
    conda:
        STAR_ENV_FILE
    input:
        STAR_GENOME_FASTA,
    output:
        directory(STAR_GENOME_DIR),
    threads: 32
    shell:
        """
        mkdir -p '{output}'
        STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir '{output}' \
        --genomeFastaFiles '{input}'
        """


rule run_star_pe_pass1:
    conda:
        STAR_ENV_FILE
    input:
        fq1_file=(
            "{sample}_R1.fastq.gz" if SKIP_TRIMMING else "{sample}_R1.trimmed.fastq.gz"
        ),
        fq1_file=(
            "{sample}_R2.fastq.gz" if SKIP_TRIMMING else "{sample}_R2.trimmed.fastq.gz"
        ),
        genome_dir=STAR_GENOME_DIR,
        gtf_file=config["ref"]["annotation"],
    output:
        STAR_PASS1_SJ_FILE,
    params:
        out_dir=join(STAR_PASS1_OUTPUT_DIR, ""),
        sjdb_overhang=lambda wildcards, input: get_sjdb_overhang(input.metadata),
    threads: config["star"]["index"]["threads"]
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --readFilesIn '{input.fq1}' '{input.fq2}' \
        --genomeDir '{input.genome_dir}' \
        --outFileNamePrefix '{params.out_dir}' \
        --outSAMattrRGline ID:{wildcards.sample} PL:Illumina SM:{wildcards.sample} LB:RNA \
        --readFilesCommand zcat \
        --sjdbGTFfile '{input.gtf_file}' \
        --sjdbOverhang {params.sjdb_overhang}
        """


rule run_star_filter_pass1_sj:
    input:
        STAR_PASS1_SJ_FILE,
    output:
        STAR_PASS1_SJ_FILTERED_FILE,
    run:
        num_filtered_novel_sj = 0
        num_novel_sj = 0
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
    conda:
        STAR_ENV_FILE
    input:
        fq1_file=(
            "{sample}_R1.fastq.gz" if SKIP_TRIMMING else "{sample}_R1.trimmed.fastq.gz"
        ),
        fq1_file=(
            "{sample}_R2.fastq.gz" if SKIP_TRIMMING else "{sample}_R2.trimmed.fastq.gz"
        ),
        genome_dir=STAR_GENOME_DIR,
        gtf_file=config["ref"]["annotation"],
        sj_file=STAR_PASS1_SJ_FILTERED_FILE,
    output:
        bam_file=STAR_PASS2_BAM_FILE,
        count_file=STAR_PASS2_READCOUNT_FILE,
    params:
        out_dir=join(STAR_PASS2_OUTPUT_DIR, ""),
        sjdb_overhang=lambda wildcards, input: get_sjdb_overhang(input.metadata),
    threads: config["star"]["align"]["threads"]
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --readFilesIn '{input.fq1}' '{input.fq2}' \
        --genomeDir '{input.genome_dir}' \
        --outFileNamePrefix '{params.out_dir}' \
        --outSAMattrRGline ID:{wildcards.sample} PL:Illumina SM:{wildcards.sample} LB:RNA \
        --readFilesCommand zcat \
        --sjdbGTFfile '{input.gtf_file}' \
        --sjdbOverhang {params.sjdb_overhang} \
        --sjdbFileChrStartEnd '{input.sj_file}' \
        """
