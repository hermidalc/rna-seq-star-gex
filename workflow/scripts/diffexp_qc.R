suppressPackageStartupMessages({
    library(Biobase)
    library(DESeq2)
    library(EDASeq)
    library(edgeR)
    library(limma)
    library(RColorBrewer)
    library(stringr)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

fig_dim <- 5
fig_res <- 300
cex_axis <- 0.4
colors <- brewer.pal(6, "Set2")
colors <- colors[1:2]

method <- snakemake@params[["method"]]
type <- snakemake@params[["type"]]
experiment <- snakemake@params[["experiment"]]
conditions <- snakemake@params[["conditions"]]
qc_legend <- snakemake@params[["qc_legend"]]
ylim <- snakemake@params[["ylim"]]

eset <- readRDS(snakemake@input[[1]])
eset <- eset[, eset$experiment == experiment]

counts <- exprs(eset)
pdata <- pData(eset)
fdata <- fData(eset)

pdata$condition <- factor(pdata$condition, levels = conditions)

if ("batch" %in% pdata && !any(is.na(pdata$batch))) {
    pdata$batch <- factor(pdata$batch)
    formula <- ~ batch + condition
} else {
    formula <- ~condition
}

design <- model.matrix(formula, data = pdata)

if (method == "counts") {
    dge <- DGEList(counts = counts, genes = fdata)
    dge <- dge[filterByExpr(dge, design), , keep.lib.sizes = FALSE]
    mat <- dge$counts
    title <- "Counts"
} else if (method %in% c("edger", "voom")) {
    dge <- DGEList(counts = counts, genes = fdata)
    dge <- dge[filterByExpr(dge, design), , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")
    # no log here as logCPMs can have negative values even with prior count
    cpms <- cpm(dge)
    mat <- cpms
    title <- "edgeR TMM"
} else if (method == "deseq2") {
    dds <- DESeqDataSetFromMatrix(counts, pdata, formula)
    dds <- dds[filterByExpr(counts, design), ]
    dds <- estimateSizeFactors(dds, quiet = TRUE)
    dds <- estimateDispersions(dds, quiet = TRUE)
    cpms <- cpm(counts(dds, normalized = TRUE))
    mat <- cpms
    title <- "DESeq2 MOR"
} else if (method == "limmarbe") {
    dge <- DGEList(counts = counts, genes = fdata)
    dge <- dge[filterByExpr(dge, design), , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")
    cpms <- cpm(dge, log = TRUE)
    cpms <- removeBatchEffect(
        cpms,
        batch = pdata$batch,
        design = model.matrix(~condition, data = pdata)
    )
    # no log as logCPMs can have negative values even with prior count
    cpms <- 2^cpms
    mat <- cpms
    title <- "limma removeBatchEffect"
}

png(
    file = snakemake@output[[1]],
    width = fig_dim, height = fig_dim, units = "in", res = fig_res
)

if (type == "rle") {
    plotRLE(
        mat,
        col = colors[as.integer(pdata$condition)],
        outline = FALSE, ylim = ylim, las = 2, cex.axis = cex_axis
    )
} else if (type == "pca") {
    plotPCA(
        mat,
        col = colors[as.integer(pdata$condition)], cex = 0.6
    )
} else if (type == "mds") {
    # plotMDS part of edgeR - pass dge directly knows what to do
    if (method == "edger") {
        mat <- dge
    }
    plotMDS(
        mat,
        col = colors[as.integer(pdata$condition)], cex = 0.6
    )
}

title(paste(experiment, title, str_to_upper(type)))
legend(
    qc_legend,
    legend = unique(pdata$condition_label),
    col = colors, pch = 15, cex = 0.8
)
invisible(dev.off())

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
