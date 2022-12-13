suppressPackageStartupMessages({
    library(Biobase)
    library(EnhancedVolcano)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

fig_dim <- 8
fig_res <- 300

method <- snakemake@params[["method"]]
experiment <- snakemake@params[["experiment"]]
conditions <- snakemake@params[["conditions"]]

fc <- snakemake@params[["fc"]]
padj <- snakemake@params[["padj"]]
padj_meth <- snakemake@params[["padj_meth"]]

lfc <- log2(fc)

eset <- readRDS(snakemake@input[[1]])
eset <- eset[, eset$experiment == experiment]

pdata <- pData(eset)

results <- read.delim(
    snakemake@input[["results"]],
    sep = "\t", header = TRUE, row.names = 1
)

if (method == "edger") {
    x <- "logFC"
    y <- "PValue"
    f <- "FDR"
    subtitle <- paste(
        "edgeR: TMM + QLFit +", ifelse(lfc > 0, "TREAT", "QLFTest")
    )
} else if (method == "deseq2") {
    x <- "log2FoldChange"
    y <- "pvalue"
    f <- "padj"
    subtitle <- paste(
        "DESeq2: MOR + nbWaldtest", ifelse(lfc > 0, "+ lfcThreshold", "")
    )
} else if (method == "voom") {
    x <- "logFC"
    y <- "P.Value"
    f <- "FDR"
    subtitle <- paste(
        "limma-voom: TMM + lmFit +", ifelse(lfc > 0, "TREAT", "eBayes")
    )
}

if ("Symbol" %in% colnames(results)) {
    labels <- results$Symbol
} else {
    labels <- row.names(results)
}

condition_labels <- unique(pdata$condition_label)
title <- paste(condition_labels[1], "vs", condition_labels[2], experiment)

png(
    file = snakemake@output[[1]],
    width = fig_dim, height = fig_dim, units = "in", res = fig_res
)
EnhancedVolcano(
    results,
    lab = labels,
    selectLab = NULL,
    x = x,
    y = y,
    xlim = c(floor(min(results[[x]])), ceiling(max(results[[x]]))),
    ylim = c(0, ceiling(max(-log10(results[[y]])))),
    pCutoff = (
        tail(results[[y]][results[[f]] < padj], n = 1)
        + head(results[[y]][results[[f]] > padj], n = 1)
    ) / 2,
    FCcutoff = lfc,
    pointSize = 3.0,
    labSize = 3.5,
    labFace = "bold",
    drawConnectors = FALSE,
    widthConnectors = 1,
    colConnectors = "grey20",
    maxoverlapsConnectors = 50,
    arrowheads = FALSE,
    boxedLabels = FALSE,
    title = title,
    subtitle = NULL,
    caption = paste("FC threshold = ", fc, " FDR = ", padj, sep = "")
)
invisible(dev.off())

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
