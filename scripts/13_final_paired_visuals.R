
library(DESeq2)
library(ggplot2)

cat("=== Final paired DEG outputs ===\n\n")

outdir <- "results/deseq2_paired/final"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

# --------------------------------------------------
# 1. Load paired analyses
# --------------------------------------------------

full <- read.delim(
    "results/deseq2_paired/paired_DESeq2_results_annotated.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

no33 <- read.delim(
    "results/deseq2_paired/paired_without33_DESeq2_results.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# --------------------------------------------------
# 2. Significant sets
# --------------------------------------------------

sig_full <- full[
    !is.na(full$padj) &
    full$padj < 0.05 &
    abs(full$log2FoldChange) >= 1,
]

sig_no33 <- no33[
    !is.na(no33$padj) &
    no33$padj < 0.05 &
    abs(no33$log2FoldChange) >= 1,
]

robust_ids <- intersect(
    sig_full$Geneid,
    sig_no33$Geneid
)

robust <- sig_full[
    sig_full$Geneid %in% robust_ids,
]

robust <- robust[
    order(robust$padj),
]

write.table(
    robust,
    file.path(
        outdir,
        "paired_robust_1066_DEGs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Paired all12 DEG:",
    nrow(sig_full), "\n")

cat("Paired without33 DEG:",
    nrow(sig_no33), "\n")

cat("Robust common DEG:",
    nrow(robust), "\n\n")

# --------------------------------------------------
# 3. Robust directions
# --------------------------------------------------

robust_up <- robust[
    robust$log2FoldChange >= 1,
]

robust_down <- robust[
    robust$log2FoldChange <= -1,
]

write.table(
    robust_up,
    file.path(
        outdir,
        "paired_robust_PD_up.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    robust_down,
    file.path(
        outdir,
        "paired_robust_PD_down.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Robust PD-up:",
    nrow(robust_up), "\n")

cat("Robust PD-down:",
    nrow(robust_down), "\n\n")

# --------------------------------------------------
# 4. Protein-coding robust set
# --------------------------------------------------

robust_pc <- robust[
    robust$gene_biotype == "protein_coding",
]

write.table(
    robust_pc,
    file.path(
        outdir,
        "paired_robust_protein_coding_DEGs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Robust protein-coding DEG:",
    nrow(robust_pc), "\n")

cat("Robust protein-coding PD-up:",
    sum(
        robust_pc$log2FoldChange >= 1
    ),
    "\n"
)

cat("Robust protein-coding PD-down:",
    sum(
        robust_pc$log2FoldChange <= -1
    ),
    "\n\n")

# --------------------------------------------------
# 5. Volcano categories
# --------------------------------------------------

full$status <- "Not significant"

full$status[
    !is.na(full$padj) &
    full$padj < 0.05 &
    full$log2FoldChange >= 1
] <- "Up in PD"

full$status[
    !is.na(full$padj) &
    full$padj < 0.05 &
    full$log2FoldChange <= -1
] <- "Down in PD"

full$minus_log10_padj <- 0

ok <- !is.na(full$padj)

full$minus_log10_padj[ok] <-
    -log10(
        pmax(
            full$padj[ok],
            .Machine$double.xmin
        )
    )

# Label top 12 robust protein-coding genes
labels <- head(
    robust_pc,
    12
)

labels$minus_log10_padj <-
    -log10(
        pmax(
            labels$padj,
            .Machine$double.xmin
        )
    )

# --------------------------------------------------
# 6. Volcano
# --------------------------------------------------

p <- ggplot(
    full,
    aes(
        x = log2FoldChange,
        y = minus_log10_padj
    )
) +
    geom_point(
        aes(color = status),
        alpha = 0.55,
        size = 1.2
    ) +
    geom_vline(
        xintercept = c(-1, 1),
        linetype = "dashed"
    ) +
    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed"
    ) +
    geom_text(
        data = labels,
        aes(
            x = log2FoldChange,
            y = minus_log10_padj,
            label = gene_name
        ),
        inherit.aes = FALSE,
        size = 3,
        vjust = -0.7
    ) +
    labs(
        title = "Paired PD vs Healthy NPC RNA-seq",
        subtitle = "DESeq2 design: ~ pair + condition",
        x = "log2 Fold Change (PD / Healthy)",
        y = "-log10 adjusted p-value",
        color = "Status"
    ) +
    theme_classic()

ggsave(
    file.path(
        outdir,
        "paired_volcano_PD_vs_Healthy.pdf"
    ),
    p,
    width = 9,
    height = 7
)

# --------------------------------------------------
# 7. Top 30 robust protein-coding heatmap
# --------------------------------------------------

vsd <- readRDS(
    "results/deseq2_paired/paired_vsd.rds"
)

mat <- assay(vsd)

top30 <- head(
    robust_pc,
    30
)

top_ids <- top30$Geneid

heat <- mat[
    top_ids,
    ,
    drop = FALSE
]

heat_z <- t(
    scale(
        t(heat)
    )
)

rownames(heat_z) <- make.unique(
    paste0(
        top30$gene_name,
        " | ",
        top30$Geneid
    )
)

meta <- as.data.frame(
    colData(vsd)
)

colnames(heat_z) <- paste(
    rownames(meta),
    paste0("Pair", meta$pair),
    meta$condition,
    sep = "_"
)

pdf(
    file.path(
        outdir,
        "paired_top30_robust_DEG_heatmap.pdf"
    ),
    width = 11,
    height = 12
)

heatmap(
    heat_z,
    scale = "none",
    margins = c(13, 14),
    cexRow = 0.7,
    cexCol = 0.65,
    main = "Top 30 robust protein-coding DEGs"
)

dev.off()

# --------------------------------------------------
# 8. Top genes
# --------------------------------------------------

cat("Top 20 robust protein-coding genes:\n\n")

print(
    head(
        robust_pc[
            ,
            c(
                "Geneid",
                "gene_name",
                "baseMean",
                "log2FoldChange",
                "padj"
            )
        ],
        20
    )
)

# --------------------------------------------------
# 9. DKK1
# --------------------------------------------------

cat("\nDKK1 final paired result:\n\n")

print(
    full[
        full$gene_name == "DKK1",
        c(
            "Geneid",
            "gene_name",
            "baseMean",
            "log2FoldChange",
            "pvalue",
            "padj"
        )
    ]
)

cat(
    "\n=== Final paired visualization completed ===\n"
)

