
OUT <- "results/ipsc_vs_npc"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. iPSC PRIMARY RESULTS
# ============================================================

ipsc <- read.delim(
    "results/deseq2_paired/paired_all_results.tsv",
    check.names = FALSE
)

# ============================================================
# 2. FIND NPC PRIMARY PAIRED RESULTS AUTOMATICALLY
# ============================================================

npc_dir <- file.path(
    Sys.getenv("HOME"),
    "bioinformatics/parkinson_npc_rnaseq/results/deseq2_paired"
)

files <- list.files(
    npc_dir,
    pattern = "\\.tsv$",
    recursive = TRUE,
    full.names = TRUE
)

cat("Searching NPC result files...\n")

candidates <- list()

for (f in files) {

    # avoid obvious non-primary / sensitivity / small summary files
    if (grepl(
        "without|sensitivity|normalized|pca|robust|heatmap|candidate|core",
        basename(f),
        ignore.case = TRUE
    )) {
        next
    }

    x <- tryCatch(
        read.delim(
            f,
            check.names = FALSE
        ),
        error = function(e) NULL
    )

    if (is.null(x)) {
        next
    }

    has_lfc <- "log2FoldChange" %in% colnames(x)

    has_gene <- any(
        c(
            "Geneid",
            "gene_id",
            "ENSEMBL",
            "ensembl"
        ) %in% colnames(x)
    )

    if (
        has_lfc &&
        has_gene &&
        nrow(x) > 10000
    ) {

        candidates[[f]] <- x
    }
}

if (length(candidates) == 0) {
    stop("Could not automatically find NPC primary paired DESeq2 results.")
}

cat("\nCandidate NPC files:\n")

for (f in names(candidates)) {
    cat(
        basename(f),
        "->",
        nrow(candidates[[f]]),
        "genes\n"
    )
}

# Prefer filename containing "all"
candidate_names <- names(candidates)

preferred <- candidate_names[
    grepl(
        "all",
        basename(candidate_names),
        ignore.case = TRUE
    )
]

if (length(preferred) > 0) {

    npc_file <- preferred[1]

} else {

    sizes <- sapply(
        candidates,
        nrow
    )

    npc_file <- names(
        which.max(sizes)
    )
}

npc <- candidates[[npc_file]]

cat("\nNPC file selected:\n")
cat(npc_file, "\n")

# ============================================================
# 3. STANDARDIZE GENE ID
# ============================================================

if ("Geneid" %in% colnames(npc)) {

    npc$gene_id <- npc$Geneid

} else if ("ENSEMBL" %in% colnames(npc)) {

    npc$gene_id <- npc$ENSEMBL

} else if ("ensembl" %in% colnames(npc)) {

    npc$gene_id <- npc$ensembl
}

if (!"gene_id" %in% colnames(npc)) {
    stop("NPC gene ID column could not be standardized.")
}

# Remove version suffix if present
ipsc$gene_id <- sub("\\..*$", "", ipsc$gene_id)
npc$gene_id  <- sub("\\..*$", "", npc$gene_id)

# ============================================================
# 4. MERGE ALL TESTED GENES
# ============================================================

ipsc_small <- data.frame(
    gene_id = ipsc$gene_id,
    iPSC_log2FC = ipsc$log2FoldChange
)

npc_small <- data.frame(
    gene_id = npc$gene_id,
    NPC_log2FC = npc$log2FoldChange
)

shared <- merge(
    ipsc_small,
    npc_small,
    by = "gene_id"
)

shared <- shared[
    complete.cases(
        shared$iPSC_log2FC,
        shared$NPC_log2FC
    ),
]

# ============================================================
# 5. GLOBAL COMPARISON
# ============================================================

global_cor <- cor(
    shared$iPSC_log2FC,
    shared$NPC_log2FC,
    method = "pearson"
)

same_direction <- (
    sign(shared$iPSC_log2FC) ==
    sign(shared$NPC_log2FC)
)

same_pct <- 100 * mean(
    same_direction
)

# Genes showing noticeable change in BOTH cell types
strong <- shared[
    abs(shared$iPSC_log2FC) >= 0.5 &
    abs(shared$NPC_log2FC) >= 0.5,
]

if (nrow(strong) > 1) {

    strong_cor <- cor(
        strong$iPSC_log2FC,
        strong$NPC_log2FC,
        method = "pearson"
    )

    strong_same <- 100 * mean(
        sign(strong$iPSC_log2FC) ==
        sign(strong$NPC_log2FC)
    )

} else {

    strong_cor <- NA
    strong_same <- NA
}

# ============================================================
# 6. SAVE
# ============================================================

write.table(
    shared,
    file.path(
        OUT,
        "all_shared_tested_genes_iPSC_NPC.tsv"
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
cat("GLOBAL iPSC vs NPC COMPARISON\n")
cat("========================================\n")

cat(
    "Shared tested genes:",
    nrow(shared),
    "\n"
)

cat(
    "Global LFC correlation:",
    round(global_cor, 4),
    "\n"
)

cat(
    "Global same-direction percentage:",
    round(same_pct, 1),
    "%\n"
)

cat(
    "Genes with |LFC| >= 0.5 in BOTH:",
    nrow(strong),
    "\n"
)

cat(
    "Strong-gene LFC correlation:",
    round(strong_cor, 4),
    "\n"
)

cat(
    "Strong genes same direction:",
    round(strong_same, 1),
    "%\n"
)

cat("========================================\n")

