
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

cat("=== FINAL PAIRED PATHWAY ANALYSIS ===\n\n")

outdir <- "results/enrichment_paired_final"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ==================================================
# 1. LOAD FINAL PAIRED RESULTS
# ==================================================

robust <- read.delim(
    "results/deseq2_paired/final/paired_robust_1066_DEGs.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

full <- read.delim(
    "results/deseq2_paired/paired_DESeq2_results_annotated.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

cat("Robust DEG:", nrow(robust), "\n")
cat("All paired tested genes:", nrow(full), "\n\n")


# ==================================================
# 2. PROTEIN-CODING SETS
# ==================================================

robust_pc <- robust[
    !is.na(robust$gene_biotype) &
    robust$gene_biotype == "protein_coding",
]

background_pc <- full[
    !is.na(full$gene_biotype) &
    full$gene_biotype == "protein_coding",
]

up_ens <- unique(
    robust_pc$Geneid[
        robust_pc$log2FoldChange >= 1
    ]
)

down_ens <- unique(
    robust_pc$Geneid[
        robust_pc$log2FoldChange <= -1
    ]
)

background_ens <- unique(
    background_pc$Geneid
)

cat("Robust protein-coding:", nrow(robust_pc), "\n")
cat("PD-up protein-coding:", length(up_ens), "\n")
cat("PD-down protein-coding:", length(down_ens), "\n")
cat("Background protein-coding:", length(background_ens), "\n\n")


# ==================================================
# 3. ENSEMBL -> ENTREZ
# ==================================================

map_ids <- function(ids) {

    suppressMessages(
        bitr(
            ids,
            fromType = "ENSEMBL",
            toType = c(
                "ENTREZID",
                "SYMBOL"
            ),
            OrgDb = org.Hs.eg.db
        )
    )
}

background_map <- map_ids(background_ens)
up_map <- map_ids(up_ens)
down_map <- map_ids(down_ens)

write.table(
    background_map,
    file.path(outdir, "background_mapping.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    up_map,
    file.path(outdir, "PD_up_mapping.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    down_map,
    file.path(outdir, "PD_down_mapping.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

background_entrez <- unique(
    as.character(background_map$ENTREZID)
)

up_entrez <- unique(
    as.character(up_map$ENTREZID)
)

down_entrez <- unique(
    as.character(down_map$ENTREZID)
)

cat("Mapped background:", length(background_entrez), "\n")
cat("Mapped PD-up:", length(up_entrez), "\n")
cat("Mapped PD-down:", length(down_entrez), "\n\n")


# ==================================================
# 4. FINAL ORA - GO BP
# ==================================================

cat("Running paired GO BP ORA - PD up...\n")

go_up <- enrichGO(
    gene = up_entrez,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
)

cat("Running paired GO BP ORA - PD down...\n")

go_down <- enrichGO(
    gene = down_entrez,
    universe = background_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
)


# ==================================================
# 5. FINAL ORA - REACTOME
# ==================================================

cat("Running paired Reactome ORA - PD up...\n")

react_up <- enrichPathway(
    gene = up_entrez,
    universe = background_entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
)

cat("Running paired Reactome ORA - PD down...\n")

react_down <- enrichPathway(
    gene = down_entrez,
    universe = background_entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    minGSSize = 10,
    maxGSSize = 500,
    readable = TRUE
)


# ==================================================
# 6. ORA OUTPUT FUNCTION
# ==================================================

save_ora <- function(obj, prefix, title) {

    df <- as.data.frame(obj)

    write.table(
        df,
        file.path(
            outdir,
            paste0(prefix, ".tsv")
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat(
        "\n",
        prefix,
        " significant terms: ",
        nrow(df),
        "\n",
        sep = ""
    )

    if (nrow(df) > 0) {

        cat("Top 10:\n")

        print(
            head(
                df[
                    ,
                    c(
                        "ID",
                        "Description",
                        "GeneRatio",
                        "BgRatio",
                        "p.adjust",
                        "Count"
                    )
                ],
                10
            )
        )

        pdf(
            file.path(
                outdir,
                paste0(
                    prefix,
                    "_dotplot.pdf"
                )
            ),
            width = 10,
            height = 8
        )

        print(
            dotplot(
                obj,
                showCategory = min(
                    20,
                    nrow(df)
                )
            ) +
            ggtitle(title)
        )

        dev.off()
    }

    cat("\n")
}


save_ora(
    go_up,
    "ORA_GO_BP_PD_up",
    "Paired ORA - GO BP - Up in PD"
)

save_ora(
    go_down,
    "ORA_GO_BP_PD_down",
    "Paired ORA - GO BP - Down in PD"
)

save_ora(
    react_up,
    "ORA_Reactome_PD_up",
    "Paired ORA - Reactome - Up in PD"
)

save_ora(
    react_down,
    "ORA_Reactome_PD_down",
    "Paired ORA - Reactome - Down in PD"
)


# ==================================================
# 7. PREPARE PAIRED GSEA RANKING
# ==================================================

rank_source <- full[
    !is.na(full$stat),
    c(
        "Geneid",
        "gene_name",
        "stat"
    )
]

rank_map <- suppressMessages(
    bitr(
        unique(rank_source$Geneid),
        fromType = "ENSEMBL",
        toType = c(
            "ENTREZID",
            "SYMBOL"
        ),
        OrgDb = org.Hs.eg.db
    )
)

rank_df <- merge(
    rank_source,
    rank_map,
    by.x = "Geneid",
    by.y = "ENSEMBL"
)

rank_df <- rank_df[
    !is.na(rank_df$ENTREZID) &
    !is.na(rank_df$stat),
]

# If one Entrez ID maps multiple times,
# retain gene with strongest absolute DESeq2 statistic

rank_df <- rank_df[
    order(
        abs(rank_df$stat),
        decreasing = TRUE
    ),
]

rank_df <- rank_df[
    !duplicated(rank_df$ENTREZID),
]

gene_list <- rank_df$stat

names(gene_list) <- as.character(
    rank_df$ENTREZID
)

gene_list <- sort(
    gene_list,
    decreasing = TRUE
)

write.table(
    rank_df,
    file.path(
        outdir,
        "paired_GSEA_ranked_genes.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "Paired GSEA ranked genes:",
    length(gene_list),
    "\n"
)

cat(
    "Duplicated Entrez IDs:",
    anyDuplicated(names(gene_list)),
    "\n\n"
)


# ==================================================
# 8. FINAL PAIRED GO GSEA
# ==================================================

set.seed(42)

cat("Running paired GO BP GSEA...\n")

gsea_go <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea"
)


# ==================================================
# 9. FINAL PAIRED REACTOME GSEA
# ==================================================

set.seed(42)

cat("Running paired Reactome GSEA...\n")

gsea_reactome <- gsePathway(
    geneList = gene_list,
    organism = "human",
    minGSSize = 10,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    seed = TRUE,
    by = "fgsea"
)


# ==================================================
# 10. GSEA OUTPUT FUNCTION
# ==================================================

save_gsea <- function(obj, prefix, title) {

    df <- as.data.frame(obj)

    sig <- df[
        !is.na(df$p.adjust) &
        df$p.adjust < 0.05,
    ]

    positive <- sig[
        sig$NES > 0,
    ]

    negative <- sig[
        sig$NES < 0,
    ]

    write.table(
        df,
        file.path(
            outdir,
            paste0(
                prefix,
                "_all.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    write.table(
        sig,
        file.path(
            outdir,
            paste0(
                prefix,
                "_significant.tsv"
            )
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat(
        "\n",
        prefix,
        "\n",
        sep = ""
    )

    cat(
        "Significant pathways:",
        nrow(sig),
        "\n"
    )

    cat(
        "Positive NES - PD enriched:",
        nrow(positive),
        "\n"
    )

    cat(
        "Negative NES - Healthy enriched:",
        nrow(negative),
        "\n\n"
    )

    if (nrow(positive) > 0) {

        positive <- positive[
            order(
                positive$NES,
                decreasing = TRUE
            ),
        ]

        cat("Top positive NES:\n")

        print(
            head(
                positive[
                    ,
                    c(
                        "ID",
                        "Description",
                        "NES",
                        "p.adjust",
                        "setSize"
                    )
                ],
                10
            )
        )
    }

    if (nrow(negative) > 0) {

        negative <- negative[
            order(
                negative$NES
            ),
        ]

        cat("\nTop negative NES:\n")

        print(
            head(
                negative[
                    ,
                    c(
                        "ID",
                        "Description",
                        "NES",
                        "p.adjust",
                        "setSize"
                    )
                ],
                10
            )
        )
    }

    if (nrow(sig) > 0) {

        pdf(
            file.path(
                outdir,
                paste0(
                    prefix,
                    "_dotplot.pdf"
                )
            ),
            width = 10,
            height = 8
        )

        print(
            dotplot(
                obj,
                showCategory = 20
            ) +
            ggtitle(title)
        )

        dev.off()
    }

    cat("\n")
}


save_gsea(
    gsea_go,
    "GSEA_GO_BP_paired",
    "Paired GSEA - GO Biological Process"
)

save_gsea(
    gsea_reactome,
    "GSEA_Reactome_paired",
    "Paired GSEA - Reactome"
)


# ==================================================
# 11. SAVE R OBJECTS
# ==================================================

saveRDS(
    go_up,
    file.path(
        outdir,
        "ORA_GO_BP_PD_up.rds"
    )
)

saveRDS(
    go_down,
    file.path(
        outdir,
        "ORA_GO_BP_PD_down.rds"
    )
)

saveRDS(
    react_up,
    file.path(
        outdir,
        "ORA_Reactome_PD_up.rds"
    )
)

saveRDS(
    react_down,
    file.path(
        outdir,
        "ORA_Reactome_PD_down.rds"
    )
)

saveRDS(
    gsea_go,
    file.path(
        outdir,
        "GSEA_GO_BP_paired.rds"
    )
)

saveRDS(
    gsea_reactome,
    file.path(
        outdir,
        "GSEA_Reactome_paired.rds"
    )
)


cat(
    "\n=== FINAL PAIRED PATHWAY ANALYSIS COMPLETED ===\n"
)

