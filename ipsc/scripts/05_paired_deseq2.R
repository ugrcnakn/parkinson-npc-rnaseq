library(DESeq2)
library(ggplot2)

# ============================================================
# 1. COUNT MATRIX
# ============================================================

count_mat <- read.delim(
    "results/counts/count_matrix.tsv",
    row.names = 1,
    check.names = FALSE
)

cat("Count matrix:\n")
print(dim(count_mat))

# ============================================================
# 2. METADATA
# ============================================================

meta <- read.delim(
    "metadata/sample_metadata.tsv",
    row.names = 1,
    check.names = FALSE
)

cat("\nMetadata:\n")
print(meta)

# ============================================================
# 3. SAMPLE MATCHING
# ============================================================

if (!all(colnames(count_mat) %in% rownames(meta))) {
    stop("ERROR: Count matrix ve metadata sample isimleri eslesmiyor.")
}

meta <- meta[colnames(count_mat), , drop = FALSE]

if (!identical(colnames(count_mat), rownames(meta))) {
    stop("ERROR: Metadata sirasi count matrix ile ayni degil.")
}

cat("\nSample order check: PASS\n")

# ============================================================
# 4. FACTORS
# Healthy reference olacak
# ============================================================

meta$condition <- factor(
    meta$condition,
    levels = c("Healthy", "PD")
)

meta$pair <- factor(meta$pair)

cat("\nCondition levels:\n")
print(levels(meta$condition))

cat("\nPair levels:\n")
print(levels(meta$pair))

# ============================================================
# 5. DESEQ2 DATASET
# ============================================================

dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_mat)),
    colData = meta,
    design = ~ pair + condition
)

# ============================================================
# 6. LOW-COUNT FILTERING
# En az 6 sample'da >=10 count
# ============================================================

genes_before <- nrow(dds)

keep <- rowSums(DESeq2::counts(dds) >= 10) >= 6

dds <- dds[keep, ]

genes_after <- nrow(dds)

cat("\nGenes before filtering:", genes_before, "\n")
cat("Genes after filtering :", genes_after, "\n")

# ============================================================
# 7. DESEQ2
# ============================================================

dds <- DESeq(dds)

cat("\nResults names:\n")
print(resultsNames(dds))

res <- results(
    dds,
    contrast = c("condition", "PD", "Healthy"),
    alpha = 0.05
)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

res_df <- res_df[
    order(res_df$padj, na.last = TRUE),
]

res_df <- res_df[
    ,
    c(
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj"
    )
]

write.table(
    res_df,
    "results/deseq2_paired/paired_all_results.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 8. SIGNIFICANT DEGs
# ============================================================

sig <- res_df[
    !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1,
]

up <- sig[
    sig$log2FoldChange > 0,
]

down <- sig[
    sig$log2FoldChange < 0,
]

write.table(
    sig,
    "results/deseq2_paired/paired_significant_DEGs.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    up,
    "results/deseq2_paired/paired_PD_up.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    down,
    "results/deseq2_paired/paired_PD_down.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 9. NORMALIZED COUNTS
# ============================================================

norm_counts <- DESeq2::counts(
    dds,
    normalized = TRUE
)

write.table(
    data.frame(
        gene_id = rownames(norm_counts),
        norm_counts
    ),
    "results/deseq2_paired/normalized_counts.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 10. PCA
# ============================================================

vsd <- vst(
    dds,
    blind = FALSE
)

pca_data <- plotPCA(
    vsd,
    intgroup = c("condition", "pair"),
    returnData = TRUE
)

percentVar <- round(
    100 * attr(pca_data, "percentVar"),
    1
)

write.table(
    pca_data,
    "results/deseq2_paired/pca_coordinates.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

p <- ggplot(
    pca_data,
    aes(
        x = PC1,
        y = PC2,
        color = condition,
        shape = pair,
        label = name
    )
) +
    geom_point(size = 4) +
    geom_text(
        vjust = -0.8,
        size = 3,
        show.legend = FALSE
    ) +
    xlab(
        paste0("PC1: ", percentVar[1], "% variance")
    ) +
    ylab(
        paste0("PC2: ", percentVar[2], "% variance")
    ) +
    ggtitle(
        "Parkinson iPSC RNA-seq PCA"
    ) +
    theme_bw()

ggsave(
    "results/deseq2_paired/PCA_paired_iPSC.pdf",
    plot = p,
    width = 9,
    height = 7
)

ggsave(
    "results/deseq2_paired/PCA_paired_iPSC.png",
    plot = p,
    width = 9,
    height = 7,
    dpi = 300
)

# ============================================================
# 11. FINAL SUMMARY
# ============================================================

cat("\n")
cat("========================================\n")
cat("FINAL PAIRED iPSC ANALYSIS SUMMARY\n")
cat("========================================\n")

cat("Samples:", ncol(count_mat), "\n")
cat("Genes before filter:", genes_before, "\n")
cat("Genes after filter :", genes_after, "\n")

cat(
    "Significant DEGs (padj<0.05 & |LFC|>=1):",
    nrow(sig),
    "\n"
)

cat("PD-up  :", nrow(up), "\n")
cat("PD-down:", nrow(down), "\n")

cat(
    "PCA variance: PC1 =",
    percentVar[1],
    "% | PC2 =",
    percentVar[2],
    "%\n"
)

cat("========================================\n")
