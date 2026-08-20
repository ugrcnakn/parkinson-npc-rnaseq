library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

OUT <- "results/enrichment"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. PRIMARY RESULTS
# ============================================================

all_res <- read.delim(
    "results/deseq2_paired/final/paired_all_results_annotated.tsv",
    check.names = FALSE
)

pair_consistency <- read.delim(
    "results/deseq2_paired/pair_consistency/primary_330_pair_consistency.tsv",
    check.names = FALSE
)

# ============================================================
# 2. ROBUST CORE
# same direction + |LFC| >= 0.5 in both twin pairs
# ============================================================

robust <- pair_consistency[
    pair_consistency$same_direction == TRUE &
    pair_consistency$both_abs_lfc_0.5 == TRUE,
]

cat("Robust core genes:", nrow(robust), "\n")

write.table(
    robust,
    file.path(OUT, "robust_263_DEGs.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 3. ADD ENTREZ IDs TO ROBUST GENES
# ============================================================

robust$entrez <- mapIds(
    org.Hs.eg.db,
    keys = robust$gene_id,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
)

# Background = ALL genes actually tested in DESeq2
background_entrez <- mapIds(
    org.Hs.eg.db,
    keys = all_res$gene_id,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
)

background_entrez <- unique(
    na.omit(as.character(background_entrez))
)

robust_up <- unique(
    na.omit(
        robust$entrez[
            robust$log2FoldChange > 0
        ]
    )
)

robust_down <- unique(
    na.omit(
        robust$entrez[
            robust$log2FoldChange < 0
        ]
    )
)

cat("Robust mapped genes:", sum(!is.na(robust$entrez)), "\n")
cat("Robust PD-up mapped :", length(robust_up), "\n")
cat("Robust PD-down mapped:", length(robust_down), "\n")
cat("Background mapped    :", length(background_entrez), "\n")

# ============================================================
# 4. GO BP ORA
# ============================================================

go_up <- enrichGO(
    gene = robust_up,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
)

go_down <- enrichGO(
    gene = robust_down,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
)

write.table(
    as.data.frame(go_up),
    file.path(OUT, "ORA_GO_BP_PD_up.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    as.data.frame(go_down),
    file.path(OUT, "ORA_GO_BP_PD_down.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 5. REACTOME ORA
# ============================================================

react_up <- enrichPathway(
    gene = robust_up,
    universe = background_entrez,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
)

react_down <- enrichPathway(
    gene = robust_down,
    universe = background_entrez,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
)

write.table(
    as.data.frame(react_up),
    file.path(OUT, "ORA_Reactome_PD_up.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    as.data.frame(react_down),
    file.path(OUT, "ORA_Reactome_PD_down.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 6. GSEA RANKING
# Use DESeq2 Wald statistic from ALL tested genes
# ============================================================

gsea_df <- all_res[
    !is.na(all_res$stat),
]

gsea_entrez <- mapIds(
    org.Hs.eg.db,
    keys = gsea_df$gene_id,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
)

gsea_df$entrez <- as.character(gsea_entrez)

gsea_df <- gsea_df[
    !is.na(gsea_df$entrez),
]

# remove duplicate Entrez IDs: retain strongest |stat|
gsea_df <- gsea_df[
    order(abs(gsea_df$stat), decreasing = TRUE),
]

gsea_df <- gsea_df[
    !duplicated(gsea_df$entrez),
]

geneList <- gsea_df$stat
names(geneList) <- gsea_df$entrez

geneList <- sort(
    geneList,
    decreasing = TRUE
)

# ============================================================
# 7. GO BP GSEA
# ============================================================

gsea_go <- gseGO(
    geneList = geneList,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
)

write.table(
    as.data.frame(gsea_go),
    file.path(OUT, "GSEA_GO_BP.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 8. REACTOME GSEA
# ============================================================

gsea_react <- gsePathway(
    geneList = geneList,
    organism = "human",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = FALSE
)

write.table(
    as.data.frame(gsea_react),
    file.path(OUT, "GSEA_Reactome.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ============================================================
# 9. DOTPLOTS
# ============================================================

save_dotplot <- function(obj, filename, title) {

    if (nrow(as.data.frame(obj)) == 0) {
        return(NULL)
    }

    p <- dotplot(
        obj,
        showCategory = 15
    ) +
        ggtitle(title)

    ggsave(
        file.path(OUT, filename),
        p,
        width = 10,
        height = 7,
        dpi = 300
    )
}

save_dotplot(
    go_up,
    "ORA_GO_BP_PD_up.png",
    "GO BP ORA - PD Up"
)

save_dotplot(
    go_down,
    "ORA_GO_BP_PD_down.png",
    "GO BP ORA - PD Down"
)

save_dotplot(
    react_up,
    "ORA_Reactome_PD_up.png",
    "Reactome ORA - PD Up"
)

save_dotplot(
    react_down,
    "ORA_Reactome_PD_down.png",
    "Reactome ORA - PD Down"
)

save_dotplot(
    gsea_go,
    "GSEA_GO_BP.png",
    "GO Biological Process GSEA"
)

save_dotplot(
    gsea_react,
    "GSEA_Reactome.png",
    "Reactome GSEA"
)

# ============================================================
# 10. SUMMARY
# ============================================================

cat("\n")
cat("========================================\n")
cat("PATHWAY ENRICHMENT SUMMARY\n")
cat("========================================\n")

cat("Robust DEGs:", nrow(robust), "\n")

cat(
    "GO ORA PD-up significant terms:",
    nrow(as.data.frame(go_up)),
    "\n"
)

cat(
    "GO ORA PD-down significant terms:",
    nrow(as.data.frame(go_down)),
    "\n"
)

cat(
    "Reactome ORA PD-up significant terms:",
    nrow(as.data.frame(react_up)),
    "\n"
)

cat(
    "Reactome ORA PD-down significant terms:",
    nrow(as.data.frame(react_down)),
    "\n"
)

cat(
    "GO GSEA significant terms:",
    nrow(as.data.frame(gsea_go)),
    "\n"
)

cat(
    "Reactome GSEA significant terms:",
    nrow(as.data.frame(gsea_react)),
    "\n"
)

cat("========================================\n")
