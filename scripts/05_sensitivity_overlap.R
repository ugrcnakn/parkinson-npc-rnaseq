
cat("=== Sensitivity DEG overlap ===\n\n")

all12 <- read.delim(
    "results/deseq2/all12_DESeq2_results.tsv",
    row.names = 1,
    check.names = FALSE
)

no33 <- read.delim(
    "results/deseq2/without_SRR16117933_DESeq2_results.tsv",
    row.names = 1,
    check.names = FALSE
)

sig_all <- rownames(
    all12[
        !is.na(all12$padj) &
        all12$padj < 0.05 &
        abs(all12$log2FoldChange) >= 1,
    ]
)

sig_no33 <- rownames(
    no33[
        !is.na(no33$padj) &
        no33$padj < 0.05 &
        abs(no33$log2FoldChange) >= 1,
    ]
)

common <- intersect(sig_all, sig_no33)
only_all <- setdiff(sig_all, sig_no33)
only_no33 <- setdiff(sig_no33, sig_all)

cat("Significant all12:", length(sig_all), "\n")
cat("Significant without33:", length(sig_no33), "\n")
cat("Common significant:", length(common), "\n")
cat("Only all12:", length(only_all), "\n")
cat("Only without33:", length(only_no33), "\n")

jaccard <- length(common) / length(union(sig_all, sig_no33))

cat(
    "Jaccard overlap:",
    round(jaccard, 4),
    "\n"
)

direction_same <- sign(all12[common, "log2FoldChange"]) ==
                  sign(no33[common, "log2FoldChange"])

cat(
    "Common DEG same direction:",
    sum(direction_same),
    "/",
    length(common),
    "\n"
)

cat(
    "Direction concordance:",
    round(mean(direction_same) * 100, 2),
    "%\n"
)

