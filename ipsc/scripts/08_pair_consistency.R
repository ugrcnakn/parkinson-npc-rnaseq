library(DESeq2)
library(ggplot2)

OUT <- "results/deseq2_paired/pair_consistency"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. INPUT
# ============================================================

counts_all <- read.delim(
    "results/counts/count_matrix.tsv",
    row.names = 1,
    check.names = FALSE
)

meta_all <- read.delim(
    "metadata/sample_metadata.tsv",
    row.names = 1,
    check.names = FALSE
)

meta_all <- meta_all[colnames(counts_all), , drop = FALSE]

meta_all$condition <- factor(
    meta_all$condition,
    levels = c("Healthy", "PD")
)

# Same global filtering rule as primary analysis
keep <- rowSums(counts_all >= 10) >= 6
counts_filt <- counts_all[keep, ]

cat("Genes used:", nrow(counts_filt), "\n")

# ============================================================
# 2. FUNCTION: ANALYZE ONE PAIR
# ============================================================

run_pair <- function(pair_number) {

    samples <- rownames(meta_all)[
        meta_all$pair == pair_number
    ]

    counts_sub <- counts_filt[, samples, drop = FALSE]
    meta_sub <- meta_all[samples, , drop = FALSE]

    meta_sub$condition <- factor(
        meta_sub$condition,
        levels = c("Healthy", "PD")
    )

    dds <- DESeqDataSetFromMatrix(
        countData = round(as.matrix(counts_sub)),
        colData = meta_sub,
        design = ~ condition
    )

    dds <- DESeq(dds, quiet = TRUE)

    res <- results(
        dds,
        contrast = c("condition", "PD", "Healthy")
    )

    df <- as.data.frame(res)

    data.frame(
        gene_id = rownames(df),
        log2FoldChange = df$log2FoldChange,
        padj = df$padj
    )
}

# ============================================================
# 3. PAIR 1 AND PAIR 2
# ============================================================

p1 <- run_pair("1")
p2 <- run_pair("2")

colnames(p1)[2:3] <- c(
    "pair1_log2FC",
    "pair1_padj"
)

colnames(p2)[2:3] <- c(
    "pair2_log2FC",
    "pair2_padj"
)

merged <- merge(
    p1,
    p2,
    by = "gene_id"
)

# ============================================================
# 4. PRIMARY 330 DEG SET
# ============================================================

primary <- read.delim(
    "results/deseq2_paired/final/paired_330_DEGs_annotated.tsv",
    check.names = FALSE
)

primary_small <- primary[
    ,
    c(
        "gene_id",
        "symbol",
        "log2FoldChange",
        "padj",
        "direction"
    )
]

combined <- merge(
    primary_small,
    merged,
    by = "gene_id"
)

combined$same_direction <- (
    sign(combined$pair1_log2FC) ==
    sign(combined$pair2_log2FC)
)

combined$both_abs_lfc_0.5 <- (
    abs(combined$pair1_log2FC) >= 0.5 &
    abs(combined$pair2_log2FC) >= 0.5
)

write.table(
    combined,
    file.path(
        OUT,
        "primary_330_pair_consistency.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

consistent <- combined[
    combined$same_direction,
]

write.table(
    consistent,
    file.path(
        OUT,
        "same_direction_DEGs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 5. CORRELATIONS
# ============================================================

valid_all <- complete.cases(
    merged$pair1_log2FC,
    merged$pair2_log2FC
)

global_cor <- cor(
    merged$pair1_log2FC[valid_all],
    merged$pair2_log2FC[valid_all],
    method = "pearson"
)

valid_deg <- complete.cases(
    combined$pair1_log2FC,
    combined$pair2_log2FC
)

deg_cor <- cor(
    combined$pair1_log2FC[valid_deg],
    combined$pair2_log2FC[valid_deg],
    method = "pearson"
)

direction_pct <- 100 * mean(
    combined$same_direction,
    na.rm = TRUE
)

strong_consistency_pct <- 100 * mean(
    combined$same_direction &
    combined$both_abs_lfc_0.5,
    na.rm = TRUE
)

# ============================================================
# 6. SCATTER PLOT
# ============================================================

p <- ggplot(
    combined,
    aes(
        x = pair1_log2FC,
        y = pair2_log2FC,
        color = direction
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dotted"
    ) +
    geom_point(
        size = 2,
        alpha = 0.7
    ) +
    theme_bw() +
    labs(
        title = "Pair-specific effect consistency",
        subtitle = "Primary 330 iPSC DEGs",
        x = "Pair 1 log2FC (PD / Healthy)",
        y = "Pair 2 log2FC (PD / Healthy)"
    )

ggsave(
    file.path(
        OUT,
        "pair1_vs_pair2_DEG_log2FC.png"
    ),
    p,
    width = 8,
    height = 7,
    dpi = 300
)

# ============================================================
# 7. SUMMARY
# ============================================================

cat("\n")
cat("========================================\n")
cat("PAIR CONSISTENCY SUMMARY\n")
cat("========================================\n")

cat(
    "Global LFC correlation:",
    round(global_cor, 4),
    "\n"
)

cat(
    "330 DEG LFC correlation:",
    round(deg_cor, 4),
    "\n"
)

cat(
    "Primary DEGs same direction:",
    sum(combined$same_direction, na.rm = TRUE),
    "/",
    nrow(combined),
    "(",
    round(direction_pct, 1),
    "% )\n"
)

cat(
    "Same direction + |LFC|>=0.5 in both pairs:",
    sum(
        combined$same_direction &
        combined$both_abs_lfc_0.5,
        na.rm = TRUE
    ),
    "/",
    nrow(combined),
    "(",
    round(strong_consistency_pct, 1),
    "% )\n"
)

cat("========================================\n")
