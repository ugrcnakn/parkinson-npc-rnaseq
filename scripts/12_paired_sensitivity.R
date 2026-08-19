
library(DESeq2)

cat("=== Paired sensitivity analysis ===\n\n")

drop_sample <- "SRR16117933"

# --------------------------------------------------
# 1. Full paired results
# --------------------------------------------------

full <- read.delim(
    "results/deseq2_paired/paired_DESeq2_results_annotated.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# --------------------------------------------------
# 2. Original counts + paired metadata
# --------------------------------------------------

counts_full <- read.delim(
    "results/counts/count_matrix.tsv",
    row.names = 1,
    check.names = FALSE
)

meta <- read.delim(
    "metadata/sample_metadata_paired.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

rownames(meta) <- meta$sample

# Use exactly the same genes as full paired model
full_dds <- readRDS(
    "results/deseq2_paired/paired_dds.rds"
)

genes_used <- rownames(full_dds)

counts_full <- counts_full[
    genes_used,
    ,
    drop = FALSE
]

# --------------------------------------------------
# 3. Remove SRR16117933
# --------------------------------------------------

keep_samples <- setdiff(
    colnames(counts_full),
    drop_sample
)

counts_sub <- counts_full[
    ,
    keep_samples,
    drop = FALSE
]

meta_sub <- meta[
    keep_samples,
    ,
    drop = FALSE
]

meta_sub$condition <- factor(
    meta_sub$condition,
    levels = c(
        "Healthy",
        "PD"
    )
)

meta_sub$pair <- factor(
    meta_sub$pair
)

stopifnot(
    identical(
        colnames(counts_sub),
        rownames(meta_sub)
    )
)

cat("Samples remaining:",
    ncol(counts_sub),
    "\n")

cat("\nDesign table:\n")

print(
    table(
        Pair = meta_sub$pair,
        Condition = meta_sub$condition
    )
)

# --------------------------------------------------
# 4. Fresh paired DESeq2 model
# --------------------------------------------------

count_matrix <- as.matrix(counts_sub)

storage.mode(
    count_matrix
) <- "integer"

dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = meta_sub,
    design = ~ pair + condition
)

dds <- DESeq(
    dds,
    minReplicatesForReplace = Inf
)

res <- results(
    dds,
    contrast = c(
        "condition",
        "PD",
        "Healthy"
    ),
    alpha = 0.05
)

res_df <- as.data.frame(res)

res_df$Geneid <- rownames(
    res_df
)

# --------------------------------------------------
# 5. Annotation
# --------------------------------------------------

anno <- read.delim(
    "results/annotation/ensembl92_gene_annotation.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

sub <- merge(
    res_df,
    anno,
    by = "Geneid",
    all.x = TRUE
)

sub <- sub[
    order(sub$padj),
]

write.table(
    sub,
    "results/deseq2_paired/paired_without33_DESeq2_results.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# --------------------------------------------------
# 6. Significant sets
# --------------------------------------------------

sig_full <- full[
    !is.na(full$padj) &
    full$padj < 0.05 &
    abs(full$log2FoldChange) >= 1,
]

sig_sub <- sub[
    !is.na(sub$padj) &
    sub$padj < 0.05 &
    abs(sub$log2FoldChange) >= 1,
]

full_ids <- sig_full$Geneid
sub_ids <- sig_sub$Geneid

common <- intersect(
    full_ids,
    sub_ids
)

only_full <- setdiff(
    full_ids,
    sub_ids
)

only_sub <- setdiff(
    sub_ids,
    full_ids
)

# --------------------------------------------------
# 7. Global LFC comparison
# --------------------------------------------------

comparison <- merge(
    full[
        ,
        c(
            "Geneid",
            "log2FoldChange"
        )
    ],
    sub[
        ,
        c(
            "Geneid",
            "log2FoldChange"
        )
    ],
    by = "Geneid",
    suffixes = c(
        "_all12",
        "_without33"
    )
)

ok <- complete.cases(
    comparison$log2FoldChange_all12,
    comparison$log2FoldChange_without33
)

lfc_cor <- cor(
    comparison$log2FoldChange_all12[ok],
    comparison$log2FoldChange_without33[ok]
)

# --------------------------------------------------
# 8. Direction concordance
# --------------------------------------------------

full_common <- sig_full[
    match(
        common,
        sig_full$Geneid
    ),
]

sub_common <- sig_sub[
    match(
        common,
        sig_sub$Geneid
    ),
]

same_direction <- sign(
    full_common$log2FoldChange
) == sign(
    sub_common$log2FoldChange
)

jaccard <- length(common) /
    length(
        union(
            full_ids,
            sub_ids
        )
    )

# --------------------------------------------------
# 9. Output
# --------------------------------------------------

cat(
    "\nPaired all12 significant:",
    length(full_ids),
    "\n"
)

cat(
    "Paired without33 significant:",
    length(sub_ids),
    "\n"
)

cat(
    "Common significant:",
    length(common),
    "\n"
)

cat(
    "Only paired all12:",
    length(only_full),
    "\n"
)

cat(
    "Only without33:",
    length(only_sub),
    "\n"
)

cat(
    "Jaccard overlap:",
    round(jaccard, 4),
    "\n"
)

cat(
    "LFC correlation:",
    round(lfc_cor, 4),
    "\n"
)

cat(
    "Direction concordance:",
    round(
        mean(same_direction) * 100,
        2
    ),
    "%\n"
)

# --------------------------------------------------
# DKK1
# --------------------------------------------------

cat("\nDKK1 sensitivity:\n")

print(
    sub[
        sub$gene_name == "DKK1",
        c(
            "Geneid",
            "gene_name",
            "baseMean",
            "log2FoldChange",
            "padj"
        )
    ]
)

cat(
    "\n=== Paired sensitivity completed ===\n"
)

