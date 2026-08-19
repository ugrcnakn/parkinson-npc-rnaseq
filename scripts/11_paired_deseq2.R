
library(DESeq2)
library(ggplot2)

cat("=== PAIRED DESeq2: PD vs Healthy ===\n\n")

dir.create(
    "results/deseq2_paired",
    recursive = TRUE,
    showWarnings = FALSE
)

# --------------------------------------------------
# Count matrix
# --------------------------------------------------

counts <- read.delim(
    "results/counts/count_matrix.tsv",
    row.names = 1,
    check.names = FALSE
)

# --------------------------------------------------
# Metadata
# --------------------------------------------------

meta <- read.delim(
    "metadata/sample_metadata_paired.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

rownames(meta) <- meta$sample

meta <- meta[
    colnames(counts),
    ,
    drop = FALSE
]

stopifnot(
    identical(
        colnames(counts),
        rownames(meta)
    )
)

meta$condition <- factor(
    meta$condition,
    levels = c("Healthy", "PD")
)

meta$pair <- factor(meta$pair)

cat("Experimental design:\n")
print(
    table(
        Pair = meta$pair,
        Condition = meta$condition
    )
)

# --------------------------------------------------
# DESeqDataSet
# --------------------------------------------------

count_matrix <- as.matrix(counts)
storage.mode(count_matrix) <- "integer"

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = meta,
    design = ~ pair + condition
)

# Same filtering strategy
keep <- rowSums(
    counts(dds) >= 10
) >= 6

dds <- dds[keep, ]

cat(
    "\nGenes after filtering:",
    nrow(dds),
    "\n"
)

# --------------------------------------------------
# DESeq2
# --------------------------------------------------

dds <- DESeq(dds)

cat("\nDesign:\n")
print(design(dds))

cat("\nResults names:\n")
print(resultsNames(dds))

res <- results(
    dds,
    contrast = c(
        "condition",
        "PD",
        "Healthy"
    ),
    alpha = 0.05
)

res <- res[
    order(res$padj),
]

res_df <- as.data.frame(res)
res_df$Geneid <- rownames(res_df)

# --------------------------------------------------
# Annotation
# --------------------------------------------------

anno <- read.delim(
    "results/annotation/ensembl92_gene_annotation.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

annotated <- merge(
    res_df,
    anno,
    by = "Geneid",
    all.x = TRUE
)

annotated <- annotated[
    order(annotated$padj),
]

write.table(
    annotated,
    "results/deseq2_paired/paired_DESeq2_results_annotated.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# --------------------------------------------------
# Significant DEG
# --------------------------------------------------

sig <- annotated[
    !is.na(annotated$padj) &
    annotated$padj < 0.05 &
    abs(annotated$log2FoldChange) >= 1,
]

write.table(
    sig,
    "results/deseq2_paired/paired_significant_DEGs.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "\nSignificant DEG:",
    nrow(sig),
    "\n"
)

cat(
    "PD-up:",
    sum(sig$log2FoldChange >= 1),
    "\n"
)

cat(
    "PD-down:",
    sum(sig$log2FoldChange <= -1),
    "\n"
)

# --------------------------------------------------
# PCA
# --------------------------------------------------

vsd <- vst(
    dds,
    blind = TRUE
)

pca <- plotPCA(
    vsd,
    intgroup = c(
        "condition",
        "pair"
    ),
    returnData = TRUE
)

percentVar <- round(
    100 * attr(
        pca,
        "percentVar"
    ),
    1
)

write.table(
    pca,
    "results/deseq2_paired/paired_PCA_data.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

pdf(
    "results/deseq2_paired/paired_PCA.pdf",
    width = 8,
    height = 6
)

print(
    plotPCA(
        vsd,
        intgroup = c(
            "condition",
            "pair"
        )
    ) +
    ggtitle(
        paste0(
            "Paired PD vs Healthy NPC PCA | PC1 ",
            percentVar[1],
            "%, PC2 ",
            percentVar[2],
            "%"
        )
    )
)

dev.off()

# --------------------------------------------------
# Compare old unpaired analysis
# --------------------------------------------------

old <- read.delim(
    "results/deseq2/all12_DESeq2_results_annotated.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

compare <- merge(
    old[
        ,
        c(
            "Geneid",
            "log2FoldChange",
            "padj"
        )
    ],
    annotated[
        ,
        c(
            "Geneid",
            "log2FoldChange",
            "padj"
        )
    ],
    by = "Geneid",
    suffixes = c(
        "_unpaired",
        "_paired"
    )
)

ok <- complete.cases(
    compare$log2FoldChange_unpaired,
    compare$log2FoldChange_paired
)

lfc_cor <- cor(
    compare$log2FoldChange_unpaired[ok],
    compare$log2FoldChange_paired[ok]
)

old_sig <- compare$Geneid[
    !is.na(compare$padj_unpaired) &
    compare$padj_unpaired < 0.05 &
    abs(compare$log2FoldChange_unpaired) >= 1
]

paired_sig <- compare$Geneid[
    !is.na(compare$padj_paired) &
    compare$padj_paired < 0.05 &
    abs(compare$log2FoldChange_paired) >= 1
]

common <- intersect(
    old_sig,
    paired_sig
)

cat(
    "\nUnpaired vs paired LFC correlation:",
    round(lfc_cor, 4),
    "\n"
)

cat(
    "Old unpaired DEG:",
    length(old_sig),
    "\n"
)

cat(
    "Paired DEG:",
    length(paired_sig),
    "\n"
)

cat(
    "Common DEG:",
    length(common),
    "\n"
)

# --------------------------------------------------
# Genes highlighted in original publication
# --------------------------------------------------

targets <- c(
    "TNF",
    "INHBA",
    "WNT7A",
    "DKK1"
)

cat(
    "\nOriginal-study candidate genes:\n\n"
)

print(
    annotated[
        annotated$gene_name %in% targets,
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

saveRDS(
    dds,
    "results/deseq2_paired/paired_dds.rds"
)

saveRDS(
    vsd,
    "results/deseq2_paired/paired_vsd.rds"
)

cat(
    "\n=== Paired analysis completed ===\n"
)

