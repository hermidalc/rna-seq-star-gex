rule diffexp_qc:
    conda:
        "../envs/diffexp_qc.yaml"
    input:
        COUNT_ESET_FILE,
    params:
        method="{de_qc_meth}",
        type="{de_qc_type}",
        experiment=lambda wc: config["diffexp"]["experiments"][
            EXPAND_PARAMS["de_exp"].index(wc.de_exp)
        ],
        contrast=config["diffexp"]["contrast"],
        contrast_label=config["diffexp"]["contrast_label"],
        has_batches=config["diffexp"]["has_batches"],
        qc_legend=config["diffexp"]["qc_legend"],
        ylim=config["diffexp"]["rle_ylim"],
    output:
        DIFFEXP_QC_PLOT_FILE,
    log:
        DIFFEXP_QC_LOG_FILE,
    script:
        "../scripts/diffexp_qc.R"


rule diffexp:
    conda:
        "../envs/diffexp.yaml"
    input:
        COUNT_ESET_FILE,
    params:
        method="{de_meth}",
        experiment=lambda wc: config["diffexp"]["experiments"][
            EXPAND_PARAMS["de_exp"].index(wc.de_exp)
        ],
        contrast=config["diffexp"]["contrast"],
        contrast_label=config["diffexp"]["contrast_label"],
        has_batches=config["diffexp"]["has_batches"],
        fc=config["diffexp"]["fc"],
        padj=config["diffexp"]["padj"],
        padj_meth=config["diffexp"]["padj_meth"],
    output:
        results=DIFFEXP_RESULTS_TABLE_FILE,
        volcano=DIFFEXP_VOLCANO_PLOT_FILE,
    log:
        DIFFEXP_RESULTS_LOG_FILE,
    script:
        "../scripts/diffexp.R"


rule heatmap:
    conda:
        "../envs/heatmap.yaml"
    input:
        COUNT_ESET_FILE,
    params:
        method="{de_meth}",
        experiment=lambda wc: config["diffexp"]["experiments"][
            EXPAND_PARAMS["de_exp"].index(wc.de_exp)
        ],
        contrast=config["diffexp"]["contrast"],
        contrast_label=config["diffexp"]["contrast_label"],
        has_batches=config["diffexp"]["has_batches"],
        fc=config["diffexp"]["fc"],
        padj=config["diffexp"]["padj"],
        padj_meth=config["diffexp"]["padj_meth"],
        fig_w=config["heatmap"]["fig_w"],
        fig_h=config["heatmap"]["fig_h"],
        seed=config["random_seed"],
    output:
        DIFFEXP_HEATMAP_PLOT_FILE,
    log:
        DIFFEXP_HEATMAP_LOG_FILE,
    script:
        "../scripts/heatmap.R"
