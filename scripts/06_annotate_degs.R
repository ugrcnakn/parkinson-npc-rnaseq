
cat("=== DEG annotation ===\n\n")

# Ensembl annotation
anno <- read.delim(
    "results/annotation/ensembl92_gene_annotation.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# Main DESeq2 analysis
res <- read.delim(
    "results/deseq2/all12_DESeq2_results.tsv",
    row.names = 1,
    check.names = FALSE
)

res$Geneid <- rownames(res)

# Merge
annotated <- merge(
    res,
    anno,
    by = "Geneid",
    all.x = TRUE
)

# Order by adjusted p-value
annotated <- annotated[
    order(annotated$padj),
]

write.table(
    annotated,
    "results/deseq2/all12_DESeq2_results_annotated.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# Significant genes
sig <- annotated[
    !is.na(annotated$padj) &
    annotated$padj < 0.05 &
    abs(annotated$log2FoldChange) >= 1,
]

write.table(
    sig,
    "results/deseq2/all12_significant_DEGs_annotated.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Total genes:", nrow(annotated), "\n")
cat("Significant DEG:", nrow(sig), "\n\n")

cat("Top 20 significant genes:\n\n")

print(
    head(
        sig[
            ,
            c(
                "Geneid",
                "gene_name",
                "gene_biotype",
                "baseMean",
                "log2FoldChange",
                "padj"
            )
        ],
        20
    )
)

