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

fig_dim <- 8
fig_res <- 300

experiment <- snakemake@params[["experiment"]]
conditions <- snakemake@params[["conditions"]]

fc <- snakemake@params[["fc"]]
padj <- snakemake@params[["padj"]]
padj_meth <- snakemake@params[["padj_meth"]]

lfc <- log2(fc)

eset <- readRDS(snakemake@input[[1]])
eset <- eset[, eset$experiment == experiment]

counts <- exprs(eset)
pdata <- pData(eset)
fdata <- fData(eset)

pdata$condition <- factor(pdata$condition, levels = conditions)

if ("batch" %in% pdata) {
    pdata$batch <- factor(pdata$batch)
    formula <- ~ batch + condition
} else {
    formula <- ~condition
}

design <- model.matrix(formula, data = pdata)

if (snakemake@params[["method"]] == "edger") {
    dge <- DGEList(counts = counts, genes = fdata)
    dge <- dge[filterByExpr(dge, design), , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")
    dge <- estimateDisp(dge, design, robust = TRUE)
    fit <- glmQLFit(dge, design, robust = TRUE)
    if (lfc > 0) {
        glt <- glmTreat(fit, coef = ncol(design), lfc = lfc)
    } else {
        glt <- glmQLFTest(fit, coef = ncol(design))
    }
    results <- as.data.frame(
        topTags(glt, n = Inf, adjust.method = padj_meth, sort.by = "PValue")
    )
} else if (snakemake@params[["method"]] == "deseq2") {
    dds <- DESeqDataSetFromMatrix(counts, pdata, formula)
    mcols(dds) <- DataFrame(mcols(dds), fdata)
    dds <- dds[filterByExpr(counts, design), ]
    dds <- DESeq(dds, quiet = TRUE)
    results <- as.data.frame(results(
        dds,
        name = tail(resultsNames(dds), n = 1), alpha = padj,
        independentFiltering = TRUE, lfcThreshold = lfc,
        altHypothesis = "greaterAbs", pAdjustMethod = padj_meth
    ))
    lfcs_results <- as.data.frame(lfcShrink(
        dds,
        coef = length(resultsNames(dds)), type = "apeglm", lfcThreshold = lfc,
        quiet = TRUE
    ))
    results$log2FoldChange <- lfcs_results$log2FoldChange
    # with DESeq2 you have to add fdata annotation columns to results manually
    # for example to add fdata column "Symbol" always add columns before
    # sorting below
    if ("Symbol" %in% colnames(fdata)) {
        results$Symbol <- mcols(dds)$Symbol
        results <- results[, c(
            "Symbol", colnames(results)[colnames(results) != "Symbol"]
        )]
    }
    if (lfc > 0 && ("svalue" %in% colnames(results))) {
        results <- results[order(results$svalue), ]
    } else {
        results <- results[order(results$pvalue), ]
        results <- results[!is.na(results$pvalue), , drop = FALSE]
        results <- results[!is.na(results$padj), , drop = FALSE]
    }
} else if (snakemake@params[["method"]] == "voom") {
    dge <- DGEList(counts = counts, genes = fdata)
    dge <- dge[filterByExpr(dge, design), , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge, method = "TMM")
    v <- voom(dge, design)
    fit <- lmFit(v, design)
    if (lfc > 0) {
        fit <- treat(fit, lfc = lfc, robust = TRUE)
        results <- topTreat(
            fit,
            coef = ncol(design), number = Inf, adjust.method = padj_meth,
            sort.by = "P"
        )
    } else {
        fit <- eBayes(fit, robust = TRUE)
        results <- topTable(
            fit,
            coef = ncol(design), number = Inf, adjust.method = padj_meth,
            sort.by = "P"
        )
    }
}

write.table(
    data.frame("ID" = row.names(results), results),
    file = snakemake@output[[1]],
    sep = "\t", quote = FALSE, row.names = FALSE
)

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
