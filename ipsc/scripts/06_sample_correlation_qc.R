library(DESeq2)
library(ggplot2)

counts <- read.delim(
    "results/counts/count_matrix.tsv",
    row.names = 1,
    check.names = FALSE
)

meta <- read.delim(
    "metadata/sample_metadata.tsv",
    row.names = 1,
    check.names = FALSE
)

meta <- meta[colnames(counts), , drop = FALSE]

meta$condition <- factor(
    meta$condition,
    levels = c("Healthy", "PD")
)

meta$pair <- factor(meta$pair)

dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = meta,
    design = ~ pair + condition
)

keep <- rowSums(DESeq2::counts(dds) >= 10) >= 6
dds <- dds[keep, ]

dds <- estimateSizeFactors(dds)

vsd <- vst(dds, blind = FALSE)

mat <- assay(vsd)

# Pearson sample correlation
cor_mat <- cor(mat, method = "pearson")

write.table(
    cor_mat,
    "results/deseq2_paired/sample_correlation_matrix.tsv",
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

# Long format for heatmap
df <- as.data.frame(as.table(cor_mat))
colnames(df) <- c("Sample1", "Sample2", "Correlation")

p <- ggplot(
    df,
    aes(
        x = Sample1,
        y = Sample2,
        fill = Correlation
    )
) +
    geom_tile() +
    geom_text(
        aes(label = sprintf("%.3f", Correlation)),
        size = 2.5
    ) +
    scale_fill_gradient(
        low = "white",
        high = "steelblue",
        limits = c(min(cor_mat), 1)
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        )
    ) +
    labs(
        title = "iPSC sample correlations",
        x = NULL,
        y = NULL
    )

ggsave(
    "results/deseq2_paired/sample_correlation_heatmap.png",
    p,
    width = 10,
    height = 9,
    dpi = 300
)

# ------------------------------------------------
# Within-cell-line replicate correlations
# ------------------------------------------------

cat("\n====================================\n")
cat("WITHIN CELL-LINE CORRELATIONS\n")
cat("====================================\n")

for (cl in unique(meta$cell_line)) {

    samples <- rownames(meta)[meta$cell_line == cl]

    sub <- cor_mat[samples, samples, drop = FALSE]

    vals <- sub[upper.tri(sub)]

    cat(
        cl,
        ": mean =",
        round(mean(vals), 4),
        "| min =",
        round(min(vals), 4),
        "| max =",
        round(max(vals), 4),
        "\n"
    )

    pairs <- which(
        upper.tri(sub),
        arr.ind = TRUE
    )

    for (i in seq_len(nrow(pairs))) {

        r <- pairs[i, 1]
        c <- pairs[i, 2]

        cat(
            "   ",
            rownames(sub)[r],
            "vs",
            colnames(sub)[c],
            "=",
            round(sub[r, c], 4),
            "\n"
        )
    }
}

cat("====================================\n")
