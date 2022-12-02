rule norm_eset:
    conda:
        "../envs/diffexp.yaml"
    input:
        COUNT_ESET_FILE,
    params:
        method="{de_meth}",
        experiment=lambda wc: config["experiments"][
            EXPAND_PARAMS["study"].index(wc.study)
        ],
        conditions=lambda wc: config["diffexp"]["conditions"][
            EXPAND_PARAMS["study"].index(wc.study)
        ],
    output:
        NORM_ESET_FILE,
    log:
        NORM_ESET_LOG,
    script:
        "../scripts/norm_eset.R"


rule diffexp_qc_plots:
    conda:
        "../envs/diffexp_qc.yaml"
    input:
        COUNT_ESET_FILE,
    params:
        method="{de_qc_meth}",
        type="{de_qc_type}",
        experiment=lambda wc: config["experiments"][
            EXPAND_PARAMS["study"].index(wc.study)
        ],
        conditions=lambda wc: config["diffexp"]["conditions"][
            EXPAND_PARAMS["study"].index(wc.study)
        ],
        qc_legend=config["diffexp"]["qc"]["legend"],
        ylim=lambda wc: config["diffexp"]["qc"]["rle_ylim"][
            EXPAND_PARAMS["study"].index(wc.study)
        ],
    output:
        DIFFEXP_QC_PLOT_FILE,
    log:
        DIFFEXP_QC_LOG_FILE,
    script:
        "../scripts/diffexp_qc.R"


rule diffexp_results:
    conda:
        "../envs/diffexp.yaml"
    input:
        COUNT_ESET_FILE,
    params:
        method="{de_meth}",
        experiment=lambda wc: config["experiments"][
            EXPAND_PARAMS["study"].index(wc.study)
        ],
        conditions=lambda wc: config["diffexp"]["conditions"][
            EXPAND_PARAMS["study"].index(wc.study)
        ],
        fc=config["diffexp"]["fc"],
        padj=config["diffexp"]["padj"],
        padj_meth=config["diffexp"]["padj_meth"],
    output:
        DIFFEXP_RESULTS_TABLE_FILE,
    log:
        DIFFEXP_RESULTS_LOG_FILE,
    script:
        "../scripts/diffexp.R"


rule diffexp_volcano:
    conda:
        "../envs/volcano.yaml"
    input:
        eset=COUNT_ESET_FILE,
        results=DIFFEXP_RESULTS_TABLE_FILE,
    params:
        method="{de_meth}",
        experiment=lambda wc: config["experiments"][
            EXPAND_PARAMS["study"].index(wc.study)
        ],
        condition=lambda wc: config["diffexp"]["conditions"][
            EXPAND_PARAMS["study"].index(wc.study)
        ],
        fc=config["diffexp"]["fc"],
        padj=config["diffexp"]["padj"],
        padj_meth=config["diffexp"]["padj_meth"],
    output:
        DIFFEXP_VOLCANO_PLOT_FILE,
    log:
        DIFFEXP_VOLCANO_LOG_FILE,
    script:
        "../scripts/volcano.R"


rule diffexp_heatmap:
    conda:
        "../envs/heatmap.yaml"
    input:
        eset=NORM_ESET_FILE,
        results=DIFFEXP_RESULTS_TABLE_FILE,
    params:
        method="{de_meth}",
        experiment=lambda wc: config["experiments"][
            EXPAND_PARAMS["study"].index(wc.study)
        ],
        conditions=lambda wc: config["diffexp"]["conditions"][
            EXPAND_PARAMS["study"].index(wc.study)
        ],
        padj=config["diffexp"]["padj"],
        seed=config["random_seed"],
        fig_w=config["diffexp"]["heatmap"]["fig_w"],
        fig_h=config["diffexp"]["heatmap"]["fig_h"],
    output:
        DIFFEXP_HEATMAP_PLOT_FILE,
    log:
        DIFFEXP_HEATMAP_LOG_FILE,
    script:
        "../scripts/heatmap.R"
