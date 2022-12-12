# rna-seq-star-gex

STAR + edgeR/DESeq2/limma-voom RNA-seq gene expression analysis Snakemake
workflow. Uses fastp for QC and trimming and GENCODE for the reference
genome and annotations. It currently supports two group experimental
conditions, but could easily be extended to support more complex experimental
designs and contrast matrices.

Snakemake rule graph:

![Snakemake rule graph](rna-seq-star-gex.svg)
