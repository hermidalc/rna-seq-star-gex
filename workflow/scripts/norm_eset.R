suppressPackageStartupMessages({
    library(Biobase)
    library(DESeq2)
    library(edgeR)
    library(limma)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

experiment <- snakemake@params[["experiment"]]
conditions <- snakemake@params[["conditions"]]

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

if (snakemake@params[["method"]] %in% c("edger", "voom")) {
    dge <- DGEList(counts = counts, genes = fdata)
    dge <- dge[filterByExpr(dge, design), , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")
    adata <- cpm(dge, log = TRUE)
} else if (snakemake@params[["method"]] == "deseq2") {
    dds <- DESeqDataSetFromMatrix(counts, pdata, formula)
    dds <- dds[filterByExpr(counts, design), ]
    dds <- estimateSizeFactors(dds, quiet = TRUE)
    dds <- estimateDispersions(dds, quiet = TRUE)
    vst <- varianceStabilizingTransformation(dds, blind = FALSE)
    adata <- assay(vst)
}

saveRDS(
    ExpressionSet(
        assayData = as.matrix(adata),
        phenoData = AnnotatedDataFrame(pdata),
        featureData = AnnotatedDataFrame(
            fdata[row.names(fdata) %in% row.names(adata), , drop = FALSE]
        ),
    ),
    snakemake@output[[1]]
)

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
