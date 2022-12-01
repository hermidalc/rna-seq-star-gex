suppressPackageStartupMessages({
    library(Biobase)
    library(ComplexHeatmap)
    library(circlize)
    library(data.table)
    library(edgeR)
    library(GetoptLong)
    library(RColorBrewer)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

fig_res <- 300
fig_w <- snakemake@params[["fig_w"]]
fig_h <- snakemake@params[["fig_h"]]

experiment <- snakemake@params[["experiment"]]
contrast <- snakemake@params[["contrast"]]
contrast_label <- snakemake@params[["contrast_label"]]
has_batches <- snakemake@params[["has_batches"]]

fc <- snakemake@params[["fc"]]
padj <- snakemake@params[["padj"]]
padj_meth <- snakemake@params[["padj_meth"]]
seed <- snakemake@params[["seed"]]

lfc <- log2(fc)
set.seed(seed)

eset <- readRDS(snakemake@input[[1]])
eset <- eset[, eset$experiment == experiment]

counts <- exprs(eset)
pdata <- pData(eset)
fdata <- fData(eset)

pdata$condition <- factor(pdata$condition, levels = contrast)

if (has_batches) {
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
    cpms <- cpm(dge, log = TRUE)
    mat <- cpms
    f <- "FDR"
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
    vst <- varianceStabilizingTransformation(dds, blind = FALSE)
    vsd <- assay(vst)
    mat <- vsd
    f <- "padj"
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
    cpms <- cpm(dge, log = TRUE)
    mat <- cpms
    f <- "FDR"
}

if (has_batches) {
    mat <- removeBatchEffect(
        mat,
        batch = pdata$batch
        # design = model.matrix(~condition, data = pdata)
    )
}

mat <- as.data.frame(mat)
if (snakemake@params[["method"]] == "deseq2") {
    mat$Symbol <- mcols(dds)$Symbol
} else {
    mat$Symbol <- dge$genes$Symbol
}
DT <- setDT(mat)
DT <- DT[DT[, .I[which.max(abs(rowSums(.SD)))], by = Symbol]$V1]
mat <- as.matrix(DT, rownames = "Symbol")
storage.mode(mat) <- "double"
mat <- mat[
    row.names(mat) %in% results$Symbol[results[[f]] < padj], ,
    drop = FALSE
]
colnames(mat) <- ifelse(
    is.na(pdata$sample_label), pdata$sample_name, pdata$sample_label
)

# z-score calc from limma::coolmap
M <- rowMeans(mat, na.rm = TRUE)
DF <- ncol(mat) - 1L
IsNA <- is.na(mat)
if (any(IsNA)) {
    mode(IsNA) <- "integer"
    DF <- DF - rowSums(IsNA)
    DF[DF == 0L] <- 1L
}
z <- mat - M
V <- rowSums(z^2L, na.rm = TRUE) / DF
z <- z / sqrt(V + 0.01)

title <- paste(contrast_label[1], "vs", contrast_label[2], experiment)

condition_colors <- c("azure3", "peachpuff")
names(condition_colors) <- contrast_label

genes <- NULL

ht_opt(
    legend_title_gp = gpar(fontsize = 9),
    legend_labels_gp = gpar(fontsize = 8),
    heatmap_column_names_gp = gpar(fontsize = 6),
    heatmap_column_title_gp = gpar(fontsize = 12),
    heatmap_row_title_gp = gpar(fontsize = 12),
    TITLE_PADDING = unit(c(4, 3), "mm")
)

ht <- Heatmap(
    z,
    name = "z-score",
    column_title = paste(title, qq("clustering of @{nrow(z)} DE genes")),
    column_names_rot = 45,
    cluster_columns = FALSE,
    clustering_distance_rows = "pearson",
    # col = colorRamp2(
    #     c(-3, 0, 3), c("skyblue4", "lightgoldenrodyellow", "red3")
    # ),
    col = colorRampPalette(rev(brewer.pal(
        n = 7, name = "RdYlBu"
    )))(100),
    heatmap_legend_param = list(
        at = c(-3, 0, 3),
        grid_height = unit(2, "mm"),
        legend_width = unit(20, "mm"),
        direction = "horizontal"
    ),
    top_annotation = HeatmapAnnotation(
        condition = factor(
            pdata$condition,
            levels = contrast, labels = contrast_label
        ),
        col = list(condition = condition_colors),
        show_annotation_name = FALSE,
        show_legend = TRUE,
        annotation_legend_param = list(nrow = 1),
        annotation_height = 1,
        height = unit(0.1, "mm"),
        annotation_name_gp = gpar(fontsize = 4)
    ),
    use_raster = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    show_heatmap_legend = TRUE,
    # row_km = 2,
    row_km_repeats = 1000
)

if (!is.null(genes)) {
    ht <- ht + rowAnnotation(
        link = anno_mark(
            at = which(row.names(z) %in% genes),
            labels = row.names(z)[row.names(z) %in% genes],
            labels_gp = gpar(fontsize = 12), padding = unit(1, "mm")
        )
    )
}

png(
    file = snakemake@output[[1]],
    width = fig_w, height = fig_h, units = "in", res = fig_res
)
draw(ht, heatmap_legend_side = "bottom", merge_legends = TRUE)
invisible(dev.off())

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
