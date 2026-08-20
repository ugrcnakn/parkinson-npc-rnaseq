
OUT <- "results/ipsc_vs_npc"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. INPUT
# ============================================================

ipsc <- read.delim(
    "results/enrichment/robust_263_DEGs.tsv",
    check.names = FALSE
)

npc <- read.delim(
    file.path(
        Sys.getenv("HOME"),
        "bioinformatics/parkinson_npc_rnaseq",
        "results/deseq2_paired/final",
        "paired_robust_1066_DEGs.tsv"
    ),
    check.names = FALSE
)

# Standardize NPC gene ID column
npc$gene_id <- npc$Geneid

# ============================================================
# 2. SHARED / UNIQUE GENES
# ============================================================

shared_ids <- intersect(
    ipsc$gene_id,
    npc$gene_id
)

ipsc_only <- setdiff(
    ipsc$gene_id,
    npc$gene_id
)

npc_only <- setdiff(
    npc$gene_id,
    ipsc$gene_id
)

# ============================================================
# 3. BUILD SHARED TABLE
# ============================================================

ipsc_small <- data.frame(
    gene_id = ipsc$gene_id,
    symbol_iPSC = ipsc$symbol,
    iPSC_log2FC = ipsc$log2FoldChange
)

npc_small <- data.frame(
    gene_id = npc$gene_id,
    symbol_NPC = npc$gene_name,
    NPC_log2FC = npc$log2FoldChange
)

shared <- merge(
    ipsc_small,
    npc_small,
    by = "gene_id"
)

# Use available symbol
shared$symbol <- ifelse(
    !is.na(shared$symbol_iPSC) &
    shared$symbol_iPSC != "",
    shared$symbol_iPSC,
    shared$symbol_NPC
)

# ============================================================
# 4. DIRECTION CONSISTENCY
# ============================================================

shared$same_direction <- (
    sign(shared$iPSC_log2FC) ==
    sign(shared$NPC_log2FC)
)

shared$direction <- ifelse(
    shared$same_direction &
    shared$iPSC_log2FC > 0,
    "PD_up_both",
    ifelse(
        shared$same_direction &
        shared$iPSC_log2FC < 0,
        "PD_down_both",
        "opposite_direction"
    )
)

# ============================================================
# 5. STATISTICS
# ============================================================

lfc_cor <- cor(
    shared$iPSC_log2FC,
    shared$NPC_log2FC,
    method = "pearson",
    use = "complete.obs"
)

same_n <- sum(
    shared$same_direction,
    na.rm = TRUE
)

same_pct <- 100 * mean(
    shared$same_direction,
    na.rm = TRUE
)

up_both <- sum(
    shared$direction == "PD_up_both"
)

down_both <- sum(
    shared$direction == "PD_down_both"
)

opposite <- sum(
    shared$direction == "opposite_direction"
)

# ============================================================
# 6. SAVE RESULTS
# ============================================================

shared <- shared[
    order(
        abs(shared$iPSC_log2FC) +
        abs(shared$NPC_log2FC),
        decreasing = TRUE
    ),
]

write.table(
    shared,
    file.path(
        OUT,
        "shared_50_robust_genes_iPSC_NPC.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    data.frame(gene_id = ipsc_only),
    file.path(
        OUT,
        "iPSC_only_213_robust_genes.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    data.frame(gene_id = npc_only),
    file.path(
        OUT,
        "NPC_only_1016_robust_genes.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 7. SUMMARY
# ============================================================

cat("\n")
cat("========================================\n")
cat("iPSC vs NPC ROBUST GENE COMPARISON\n")
cat("========================================\n")

cat("iPSC robust genes :", length(unique(ipsc$gene_id)), "\n")
cat("NPC robust genes  :", length(unique(npc$gene_id)), "\n")
cat("Shared genes      :", nrow(shared), "\n")
cat("iPSC-only genes   :", length(ipsc_only), "\n")
cat("NPC-only genes    :", length(npc_only), "\n")

cat(
    "Shared LFC correlation:",
    round(lfc_cor, 4),
    "\n"
)

cat(
    "Shared genes same direction:",
    same_n,
    "/",
    nrow(shared),
    "(",
    round(same_pct, 1),
    "% )\n"
)

cat("PD-up in both     :", up_both, "\n")
cat("PD-down in both   :", down_both, "\n")
cat("Opposite direction:", opposite, "\n")

cat("\nTop shared genes:\n")

print(
    head(
        shared[
            ,
            c(
                "gene_id",
                "symbol",
                "iPSC_log2FC",
                "NPC_log2FC",
                "direction"
            )
        ],
        25
    ),
    row.names = FALSE
)

cat("========================================\n")

