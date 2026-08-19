
library(DESeq2)

cat("=== Differential Expression: PD vs Healthy ===\n\n")

# --------------------------------------------------
# Function for DESeq2 analysis
# --------------------------------------------------

run_analysis <- function(dds, label) {

    cat("\n===============================\n")
    cat("Analysis:", label, "\n")
    cat("===============================\n")

    dds <- DESeq(
        dds,
        minReplicatesForReplace = Inf
    )

    res <- results(
        dds,
        contrast = c("condition", "PD", "Healthy"),
        alpha = 0.05
    )

    res <- res[order(res$padj), ]

    outfile <- paste0(
        "results/deseq2/",
        label,
        "_DESeq2_results.tsv"
    )

    write.table(
        as.data.frame(res),
        outfile,
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    sig <- subset(
        as.data.frame(res),
        !is.na(padj) &
        padj < 0.05 &
        abs(log2FoldChange) >= 1
    )

    sigfile <- paste0(
        "results/deseq2/",
        label,
        "_significant_DEGs.tsv"
    )

    write.table(
        sig,
        sigfile,
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )

    cat("Tested genes:", nrow(res), "\n")
    cat(
        "padj < 0.05 & |log2FC| >= 1:",
        nrow(sig),
        "\n"
    )

    cat(
        "Upregulated in PD:",
        sum(sig$log2FoldChange > 1),
        "\n"
    )

    cat(
        "Downregulated in PD:",
        sum(sig$log2FoldChange < -1),
        "\n"
    )

    cat("\nTop 10 genes:\n")

    print(
        head(
            as.data.frame(res)[
                ,
                c(
                    "baseMean",
                    "log2FoldChange",
                    "pvalue",
                    "padj"
                )
            ],
            10
        )
    )

    return(res)
}


# --------------------------------------------------
# Load filtered dataset
# --------------------------------------------------

dds <- readRDS(
    "results/deseq2/dds_prefiltered_normalized.rds"
)

dds$condition <- relevel(
    dds$condition,
    ref = "Healthy"
)


# --------------------------------------------------
# Analysis A: all 12 samples
# --------------------------------------------------

res_all <- run_analysis(
    dds,
    "all12"
)


# --------------------------------------------------
# Analysis B: remove SRR16117933
# --------------------------------------------------

dds_no33 <- dds[
    ,
    colnames(dds) != "SRR16117933"
]

dds_no33$condition <- droplevels(
    dds_no33$condition
)

res_no33 <- run_analysis(
    dds_no33,
    "without_SRR16117933"
)


# --------------------------------------------------
# Compare log2 fold changes
# --------------------------------------------------

common <- intersect(
    rownames(res_all),
    rownames(res_no33)
)

comparison <- data.frame(

    Geneid = common,

    LFC_all12 =
        res_all[common, "log2FoldChange"],

    padj_all12 =
        res_all[common, "padj"],

    LFC_without33 =
        res_no33[common, "log2FoldChange"],

    padj_without33 =
        res_no33[common, "padj"]
)

comparison$LFC_difference <-
    comparison$LFC_without33 -
    comparison$LFC_all12

write.table(
    comparison,
    "results/deseq2/all12_vs_without33.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# --------------------------------------------------
# Global LFC correlation
# --------------------------------------------------

ok <- complete.cases(
    comparison$LFC_all12,
    comparison$LFC_without33
)

lfc_cor <- cor(
    comparison$LFC_all12[ok],
    comparison$LFC_without33[ok],
    method = "pearson"
)

cat("\n===============================\n")
cat("Sensitivity comparison\n")
cat("===============================\n")

cat(
    "LFC correlation all12 vs without33:",
    round(lfc_cor, 4),
    "\n"
)

cat("\nAnalysis completed.\n")

