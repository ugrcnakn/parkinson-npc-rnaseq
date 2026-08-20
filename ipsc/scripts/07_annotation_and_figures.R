library(org.Hs.eg.db)
library(AnnotationDbi)
library(DESeq2)
library(ggplot2)

OUT <- "results/deseq2_paired/final"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. DESeq2 results
# ============================================================

res <- read.delim(
    "results/deseq2_paired/paired_all_results.tsv",
    check.names = FALSE
)

# Ensembl -> gene annotation
res$symbol <- mapIds(
    org.Hs.eg.db,
    keys = res$gene_id,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
)

res$entrez <- mapIds(
    org.Hs.eg.db,
    keys = res$gene_id,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
)

res$gene_name <- mapIds(
    org.Hs.eg.db,
    keys = res$gene_id,
    column = "GENENAME",
    keytype = "ENSEMBL",
    multiVals = "first"
)

res$direction <- "Not_significant"

res$direction[
    !is.na(res$padj) &
    res$padj < 0.05 &
    res$log2FoldChange >= 1
] <- "PD_up"

res$direction[
    !is.na(res$padj) &
    res$padj < 0.05 &
    res$log2FoldChange <= -1
] <- "PD_down"

sig <- res[res$direction != "Not_significant", ]

write.table(
    res,
    file.path(OUT, "paired_all_results_annotated.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    sig,
    file.path(OUT, "paired_330_DEGs_annotated.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 2. VOLCANO
# ============================================================

res$minus_log10_padj <- -log10(res$padj)

res$label <- ""

top_label <- sig[
    order(sig$padj),
]

top_label <- head(top_label, 15)

res$label[
    res$gene_id %in% top_label$gene_id
] <- ifelse(
    is.na(res$symbol[
        res$gene_id %in% top_label$gene_id
    ]),
    res$gene_id[
        res$gene_id %in% top_label$gene_id
    ],
    res$symbol[
        res$gene_id %in% top_label$gene_id
    ]
)

p <- ggplot(
    res,
    aes(
        x = log2FoldChange,
        y = minus_log10_padj
    )
) +
    geom_point(
        aes(color = direction),
        alpha = 0.65,
        size = 1.4
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
        data = res[res$label != "", ],
        aes(label = label),
        size = 3,
        vjust = -0.5,
        show.legend = FALSE
    ) +
    scale_color_manual(
        values = c(
            "PD_down" = "#377EB8",
            "Not_significant" = "grey70",
            "PD_up" = "#E41A1C"
        )
    ) +
    theme_bw() +
    labs(
        title = "iPSC PD vs Healthy - Paired DESeq2",
        subtitle = "padj < 0.05 and |log2FC| >= 1",
        x = "log2 Fold Change (PD / Healthy)",
        y = "-log10 adjusted p-value",
        color = NULL
    )

ggsave(
    file.path(OUT, "paired_volcano_iPSC.png"),
    p,
    width = 9,
    height = 7,
    dpi = 300
)

ggsave(
    file.path(OUT, "paired_volcano_iPSC.pdf"),
    p,
    width = 9,
    height = 7
)

# ============================================================
# 3. TOP 30 HEATMAP
# ============================================================

count_mat <- read.delim(
    "results/counts/count_matrix.tsv",
    row.names = 1,
    check.names = FALSE
)

meta <- read.delim(
    "metadata/sample_metadata.tsv",
    row.names = 1,
    check.names = FALSE
)

meta <- meta[colnames(count_mat), , drop = FALSE]

meta$condition <- factor(
    meta$condition,
    levels = c("Healthy", "PD")
)

meta$pair <- factor(meta$pair)

dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(count_mat)),
    colData = meta,
    design = ~ pair + condition
)

keep <- rowSums(counts(dds) >= 10) >= 6
dds <- dds[keep, ]

dds <- estimateSizeFactors(dds)

vsd <- vst(dds, blind = FALSE)

top30 <- head(sig[order(sig$padj), ], 30)

mat <- assay(vsd)[top30$gene_id, , drop = FALSE]

# row z-score
mat_z <- t(scale(t(mat)))

labels <- top30$symbol

labels[
    is.na(labels) | labels == ""
] <- top30$gene_id[
    is.na(labels) | labels == ""
]

rownames(mat_z) <- make.unique(labels)

df <- as.data.frame(as.table(mat_z))
colnames(df) <- c("Gene", "Sample", "Zscore")

df$Sample <- factor(
    df$Sample,
    levels = colnames(mat_z)
)

p_heat <- ggplot(
    df,
    aes(
        x = Sample,
        y = Gene,
        fill = Zscore
    )
) +
    geom_tile() +
    scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),
        axis.text.y = element_text(
            size = 7
        )
    ) +
    labs(
        title = "Top 30 iPSC Differentially Expressed Genes",
        x = NULL,
        y = NULL,
        fill = "Z-score"
    )

ggsave(
    file.path(OUT, "paired_top30_DEG_heatmap.png"),
    p_heat,
    width = 10,
    height = 10,
    dpi = 300
)

ggsave(
    file.path(OUT, "paired_top30_DEG_heatmap.pdf"),
    p_heat,
    width = 10,
    height = 10
)

# ============================================================
# 4. SUMMARY
# ============================================================

cat("\n====================================\n")
cat("ANNOTATION SUMMARY\n")
cat("====================================\n")

cat("Total tested genes:", nrow(res), "\n")
cat("Significant DEGs :", nrow(sig), "\n")
cat(
    "Annotated with SYMBOL:",
    sum(!is.na(sig$symbol)),
    "\n"
)
cat(
    "PD-up:",
    sum(sig$direction == "PD_up"),
    "\n"
)
cat(
    "PD-down:",
    sum(sig$direction == "PD_down"),
    "\n"
)

cat("\nTop 15 annotated DEGs:\n")

print(
    top_label[
        ,
        c(
            "gene_id",
            "symbol",
            "log2FoldChange",
            "padj"
        )
    ],
    row.names = FALSE
)

cat("====================================\n")
