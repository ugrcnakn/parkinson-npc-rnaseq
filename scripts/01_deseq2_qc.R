
library(DESeq2)

cat("=== Parkinson NPC RNA-seq: DESeq2 QC ===\n\n")

# --------------------------------------------------
# 1. Count matrix
# --------------------------------------------------

count_data <- read.delim(
    "results/counts/count_matrix.tsv",
    row.names = 1,
    check.names = FALSE
)

cat("Count matrix:\n")
cat("Genes:", nrow(count_data), "\n")
cat("Samples:", ncol(count_data), "\n\n")


# --------------------------------------------------
# 2. Sample metadata
# --------------------------------------------------

meta <- read.delim(
    "metadata/sample_metadata.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

rownames(meta) <- meta$sample

cat("Metadata samples:", nrow(meta), "\n\n")


# --------------------------------------------------
# 3. Sample matching
# --------------------------------------------------

if (!setequal(colnames(count_data), rownames(meta))) {
    stop("ERROR: Count matrix ve metadata sample isimleri eşleşmiyor!")
}

meta <- meta[colnames(count_data), , drop = FALSE]

if (!identical(colnames(count_data), rownames(meta))) {
    stop("ERROR: Sample sırası eşleşmiyor!")
}

cat("Sample matching: OK\n\n")


# --------------------------------------------------
# 4. Factors
# --------------------------------------------------

meta$condition <- factor(
    meta$condition,
    levels = c("Healthy", "PD")
)

meta$line_code <- factor(meta$line_code)

cat("Condition dağılımı:\n")
print(table(meta$condition))

cat("\nLine code dağılımı:\n")
print(table(meta$line_code))


# --------------------------------------------------
# 5. Integer count matrix
# --------------------------------------------------

count_matrix <- as.matrix(count_data)
storage.mode(count_matrix) <- "integer"


# --------------------------------------------------
# 6. DESeqDataSet
# --------------------------------------------------

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData   = meta,
    design    = ~ condition
)

cat("\nGenes before filtering:", nrow(dds), "\n")


# --------------------------------------------------
# 7. Low-count filtering
#
# Keep genes with >=10 reads in at least 6 samples
# --------------------------------------------------

keep <- rowSums(counts(dds) >= 10) >= 6

dds <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")
cat("Genes removed:", sum(!keep), "\n\n")


# --------------------------------------------------
# 8. Size-factor normalization
# --------------------------------------------------

dds <- estimateSizeFactors(dds)

cat("Size factors:\n")
print(sizeFactors(dds))


# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

write.table(
    normalized_counts,
    file = "results/deseq2/normalized_counts.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = NA
)


# --------------------------------------------------
# 9. Variance Stabilizing Transformation
# --------------------------------------------------

vsd <- vst(dds, blind = TRUE)


# --------------------------------------------------
# 10. PCA
# --------------------------------------------------

pca_data <- plotPCA(
    vsd,
    intgroup = c("condition", "line_code"),
    returnData = TRUE
)

percent_var <- round(
    100 * attr(pca_data, "percentVar"),
    1
)

write.table(
    pca_data,
    file = "results/deseq2/pca_data.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\nPCA variance:\n")
cat("PC1:", percent_var[1], "%\n")
cat("PC2:", percent_var[2], "%\n\n")

cat("PCA coordinates:\n")
print(
    pca_data[
        ,
        c("name", "condition", "line_code", "PC1", "PC2")
    ]
)


# --------------------------------------------------
# 11. PCA PDF
# --------------------------------------------------

pdf(
    "results/deseq2/PCA_condition_line.pdf",
    width = 8,
    height = 6
)

print(
    plotPCA(
        vsd,
        intgroup = c("condition", "line_code")
    ) +
    ggplot2::ggtitle(
        paste0(
            "Parkinson NPC RNA-seq PCA\n",
            "PC1: ", percent_var[1], "%   ",
            "PC2: ", percent_var[2], "%"
        )
    )
)

dev.off()


# --------------------------------------------------
# 12. Save objects
# --------------------------------------------------

saveRDS(
    dds,
    "results/deseq2/dds_prefiltered_normalized.rds"
)

saveRDS(
    vsd,
    "results/deseq2/vsd.rds"
)

cat("\nDESeq2 QC completed successfully.\n")

