
library(DESeq2)
library(ggplot2)

cat("=== DEG visualization and robust set ===\n\n")

# --------------------------------------------------
# 1. Load annotation and results
# --------------------------------------------------

anno <- read.delim(
    "results/annotation/ensembl92_gene_annotation.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

all12 <- read.delim(
    "results/deseq2/all12_DESeq2_results_annotated.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

no33 <- read.delim(
    "results/deseq2/without_SRR16117933_DESeq2_results.tsv",
    row.names = 1,
    check.names = FALSE
)

# --------------------------------------------------
# 2. Significant genes
# --------------------------------------------------

sig_all <- all12[
    !is.na(all12$padj) &
    all12$padj < 0.05 &
    abs(all12$log2FoldChange) >= 1,
]

sig_no33_ids <- rownames(
    no33[
        !is.na(no33$padj) &
        no33$padj < 0.05 &
        abs(no33$log2FoldChange) >= 1,
    ]
)

# --------------------------------------------------
# 3. Robust DEG set
# --------------------------------------------------

robust <- sig_all[
    sig_all$Geneid %in% sig_no33_ids,
]

robust <- robust[
    order(robust$padj),
]

write.table(
    robust,
    "results/deseq2/robust_1099_DEGs_annotated.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Main significant DEG:", nrow(sig_all), "\n")
cat("Robust common DEG:", nrow(robust), "\n\n")

# --------------------------------------------------
# 4. Protein-coding DEGs
# --------------------------------------------------

protein_sig <- sig_all[
    sig_all$gene_biotype == "protein_coding",
]

protein_sig <- protein_sig[
    order(protein_sig$padj),
]

write.table(
    protein_sig,
    "results/deseq2/all12_protein_coding_DEGs.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Protein-coding significant DEG:",
    nrow(protein_sig), "\n")

cat(
    "Protein-coding PD-up:",
    sum(protein_sig$log2FoldChange > 1),
    "\n"
)

cat(
    "Protein-coding PD-down:",
    sum(protein_sig$log2FoldChange < -1),
    "\n\n"
)

# --------------------------------------------------
# 5. Volcano categories
# --------------------------------------------------

all12$status <- "Not significant"

all12$status[
    !is.na(all12$padj) &
    all12$padj < 0.05 &
    all12$log2FoldChange >= 1
] <- "Up in PD"

all12$status[
    !is.na(all12$padj) &
    all12$padj < 0.05 &
    all12$log2FoldChange <= -1
] <- "Down in PD"

all12$minus_log10_padj <- 0

ok <- !is.na(all12$padj)

all12$minus_log10_padj[ok] <-
    -log10(
        pmax(
            all12$padj[ok],
            .Machine$double.xmin
        )
    )

# Top protein-coding genes for labels
label_data <- head(
    protein_sig,
    12
)

label_data$minus_log10_padj <-
    -log10(
        pmax(
            label_data$padj,
            .Machine$double.xmin
        )
    )

# --------------------------------------------------
# 6. Volcano plot
# --------------------------------------------------

p <- ggplot(
    all12,
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
        data = label_data,
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
        title = "PD vs Healthy — NPC RNA-seq",
        subtitle = "DESeq2: padj < 0.05 and |log2FC| >= 1",
        x = "log2 Fold Change (PD / Healthy)",
        y = "-log10 adjusted p-value",
        color = "Status"
    ) +
    theme_classic()

ggsave(
    "results/deseq2/volcano_PD_vs_Healthy.pdf",
    p,
    width = 9,
    height = 7
)

# --------------------------------------------------
# 7. Top 30 protein-coding DEG heatmap
# --------------------------------------------------

vsd <- readRDS(
    "results/deseq2/vsd.rds"
)

mat <- assay(vsd)

top30 <- head(
    protein_sig,
    30
)

top_ids <- top30$Geneid

heat <- mat[
    top_ids,
    ,
    drop = FALSE
]

# Z-score each gene
heat_z <- t(
    scale(
        t(heat)
    )
)

gene_labels <- paste0(
    top30$gene_name,
    " | ",
    top30$Geneid
)

rownames(heat_z) <- make.unique(gene_labels)

meta <- as.data.frame(
    colData(vsd)
)

sample_labels <- paste(
    rownames(meta),
    meta$condition,
    meta$line_code,
    sep = "_"
)

colnames(heat_z) <- sample_labels

pdf(
    "results/deseq2/top30_DEG_heatmap.pdf",
    width = 11,
    height = 12
)

heatmap(
    heat_z,
    scale = "none",
    margins = c(13, 14),
    cexRow = 0.7,
    cexCol = 0.7,
    main = "Top 30 protein-coding DEGs"
)

dev.off()

# --------------------------------------------------
# 8. Print top genes
# --------------------------------------------------

cat("\nTop 15 protein-coding DEG:\n\n")

print(
    head(
        protein_sig[
            ,
            c(
                "Geneid",
                "gene_name",
                "baseMean",
                "log2FoldChange",
                "padj"
            )
        ],
        15
    )
)

cat("\nVisualization completed.\n")

