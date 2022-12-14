suppressPackageStartupMessages({
    library(data.table)
    library(fgsea)
    library(msigdbr)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

fig_dim <- 8
fig_res <- 300

de_meth <- snakemake@params[["de_meth"]]
species <- snakemake@params[["species"]]
seed <- snakemake@params[["seed"]]
gsea_padj <- snakemake@params[["gsea_padj"]]
db_cat <- snakemake@params[["db_cat"]]
db_sub <- snakemake@params[["db_sub"]]
collapse <- snakemake@params[["collapse"]]

results <- read.delim(
    snakemake@input[[1]],
    sep = "\t", header = TRUE, row.names = 1
)

if (de_meth == "edger") {
    lfc <- "logFC"
    pval <- "PValue"
} else if (de_meth == "deseq2") {
    lfc <- "log2FoldChange"
    pval <- "pvalue"
} else if (de_meth == "voom") {
    lfc <- "logFC"
    pval <- "P.Value"
}

ranks <- results[[lfc]] * -log10(results[[pval]])

msigdb_ranks <- ranks
if ("Symbol" %in% colnames(results)) {
    names(msigdb_ranks) <- results$Symbol
} else {
    names(msigdb_ranks) <- row.names(results)
}
msigdb_ranks <- msigdb_ranks[order(abs(msigdb_ranks), decreasing = TRUE)]
msigdb_ranks <- msigdb_ranks[!duplicated(names(msigdb_ranks))]
msigdb_ranks <- msigdb_ranks[order(msigdb_ranks, decreasing = TRUE)]

msigdbr_df <- msigdbr(
    species = species, category = db_cat, subcategory = db_sub
)
msigdb_pathways <- split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

set.seed(seed)
fgsea_res <- fgsea(
    pathways = msigdb_pathways, stats = msigdb_ranks,
    minSize = 15, maxSize = 500, eps = 0.0, nPermSimple = 10000, nproc = 1
)
fgsea_res <- fgsea_res[order(pval)][padj < gsea_padj]

if (collapse) {
    msigdb_collapsed_pathways <- collapsePathways(
        fgsea_res, msigdb_pathways, msigdb_ranks
    )
    fgsea_res <- fgsea_res[
        pathway %in% msigdb_collapsed_pathways$mainPathways
    ]
}

fwrite(
    fgsea_res,
    file = snakemake@output[[1]],
    sep = "\t", sep2 = c("", " ", "")
)

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
