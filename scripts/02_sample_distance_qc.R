
library(DESeq2)

cat("=== Sample distance / correlation QC ===\n\n")

vsd <- readRDS("results/deseq2/vsd.rds")

mat <- assay(vsd)
meta <- as.data.frame(colData(vsd))

# Sample-to-sample distances
sample_dist <- dist(t(mat))
dist_matrix <- as.matrix(sample_dist)

write.table(
    dist_matrix,
    "results/deseq2/sample_distance_matrix.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

# Sample correlations
cor_matrix <- cor(mat, method = "pearson")

write.table(
    cor_matrix,
    "results/deseq2/sample_correlation_matrix.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

# Labels
labels <- paste(
    rownames(meta),
    meta$condition,
    meta$line_code,
    sep = "_"
)

# Distance heatmap
pdf(
    "results/deseq2/sample_distance_heatmap.pdf",
    width = 10,
    height = 9
)

heatmap(
    dist_matrix,
    symm = TRUE,
    labRow = labels,
    labCol = labels,
    margins = c(12,12),
    main = "Sample-to-sample distance"
)

dev.off()

# Correlation heatmap
pdf(
    "results/deseq2/sample_correlation_heatmap.pdf",
    width = 10,
    height = 9
)

heatmap(
    cor_matrix,
    symm = TRUE,
    labRow = labels,
    labCol = labels,
    margins = c(12,12),
    main = "Sample correlation"
)

dev.off()

# Nearest neighbour
cat("Nearest neighbour for each sample:\n\n")

for (i in seq_len(nrow(dist_matrix))) {

    x <- dist_matrix[i,]
    x[i] <- Inf

    nearest <- which.min(x)

    cat(
        rownames(dist_matrix)[i],
        " -> ",
        colnames(dist_matrix)[nearest],
        "   distance = ",
        round(x[nearest], 2),
        "\n",
        sep = ""
    )
}

# Within-line correlations
cat("\nMean within-line correlations:\n\n")

for (LINE in levels(meta$line_code)) {

    samples <- rownames(meta)[meta$line_code == LINE]

    x <- cor_matrix[samples, samples, drop = FALSE]

    values <- x[upper.tri(x)]

    cat(
        LINE,
        ": ",
        round(mean(values), 4),
        "\n",
        sep = ""
    )
}

cat("\nQC completed.\n")

