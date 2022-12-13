suppressPackageStartupMessages({
    library(Biobase)
    library(ComplexHeatmap)
    library(circlize)
    library(data.table)
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

method <- snakemake@params[["method"]]
experiment <- snakemake@params[["experiment"]]
conditions <- snakemake@params[["conditions"]]

padj <- snakemake@params[["padj"]]
seed <- snakemake@params[["seed"]]

eset <- readRDS(snakemake@input[["eset"]])

adata <- exprs(eset)
pdata <- pData(eset)
fdata <- fData(eset)

results <- read.delim(
    snakemake@input[["results"]],
    sep = "\t", header = TRUE, row.names = 1
)

if (method %in% c("edger", "voom")) {
    f <- "FDR"
} else if (method == "deseq2") {
    f <- "padj"
}

if ("batch" %in% pdata && !any(is.na(pdata$batch))) {
    adata <- removeBatchEffect(
        adata,
        batch = pdata$batch
        # design = model.matrix(~condition, data = pdata)
    )
}

adata <- as.data.frame(adata)
adata$Symbol <- fdata$Symbol
DT <- setDT(adata)
DT <- DT[DT[, .I[which.max(abs(rowSums(.SD)))], by = Symbol]$V1]
adata <- as.matrix(DT, rownames = "Symbol")
storage.mode(adata) <- "double"
adata <- adata[
    row.names(adata) %in% results$Symbol[results[[f]] < padj], ,
    drop = FALSE
]
colnames(adata) <- ifelse(
    is.na(pdata$sample_label), pdata$sample_name, pdata$sample_label
)

# z-score calc from limma::coolmap
M <- rowMeans(adata, na.rm = TRUE)
DF <- ncol(adata) - 1L
IsNA <- is.na(adata)
if (any(IsNA)) {
    mode(IsNA) <- "integer"
    DF <- DF - rowSums(IsNA)
    DF[DF == 0L] <- 1L
}
z <- adata - M
V <- rowSums(z^2L, na.rm = TRUE) / DF
z <- z / sqrt(V + 0.01)

condition_labels <- unique(pdata$condition_label)
title <- paste(condition_labels[1], "vs", condition_labels[2], experiment)

condition_colors <- c("azure3", "peachpuff")
names(condition_colors) <- condition_labels

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
            levels = conditions, labels = condition_labels
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
