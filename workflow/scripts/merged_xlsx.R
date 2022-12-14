suppressPackageStartupMessages({
    library(openxlsx)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

xlsx_data <- lapply(
    snakemake@input,
    function(f) read.delim(f, sep = "\t", header = TRUE)
)
names(xlsx_data) <- snakemake@params[["names"]]

write.xlsx(xlsx_data, snakemake@output[[1]])

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
