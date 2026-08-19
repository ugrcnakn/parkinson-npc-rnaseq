
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

cat("=== GSEA validation: PD vs Healthy ===\n\n")

dir.create(
    "results/gsea",
    recursive = TRUE,
    showWarnings = FALSE
)

# --------------------------------------------------
# 1. Load ALL tested DESeq2 genes
# --------------------------------------------------

res <- read.delim(
    "results/deseq2/all12_DESeq2_results_annotated.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

res <- res[
    !is.na(res$stat),
]

cat("Genes with DESeq2 statistic:",
    nrow(res), "\n")


# --------------------------------------------------
# 2. ENSEMBL -> ENTREZ
# --------------------------------------------------

id_map <- suppressMessages(
    bitr(
        unique(res$Geneid),
        fromType = "ENSEMBL",
        toType = c(
            "ENTREZID",
            "SYMBOL"
        ),
        OrgDb = org.Hs.eg.db
    )
)

rank_df <- merge(
    res[
        ,
        c(
            "Geneid",
            "gene_name",
            "stat"
        )
    ],
    id_map,
    by.x = "Geneid",
    by.y = "ENSEMBL"
)

rank_df <- rank_df[
    !is.na(rank_df$ENTREZID) &
    !is.na(rank_df$stat),
]


# --------------------------------------------------
# 3. Remove duplicate Entrez IDs
#
# Keep mapping with strongest absolute DESeq2 statistic
# --------------------------------------------------

rank_df <- rank_df[
    order(
        abs(rank_df$stat),
        decreasing = TRUE
    ),
]

rank_df <- rank_df[
    !duplicated(rank_df$ENTREZID),
]


# --------------------------------------------------
# 4. Ranked gene list
# --------------------------------------------------

gene_list <- rank_df$stat

names(gene_list) <-
    as.character(rank_df$ENTREZID)

gene_list <- sort(
    gene_list,
    decreasing = TRUE
)

cat("Unique ranked Entrez genes:",
    length(gene_list), "\n")

cat("Duplicated Entrez IDs:",
    anyDuplicated(names(gene_list)), "\n\n")

write.table(
    rank_df,
    "results/gsea/GSEA_ranked_genes.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# --------------------------------------------------
# 5. GO Biological Process GSEA
# --------------------------------------------------

set.seed(42)

cat("Running GO BP GSEA...\n")

gsea_go <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    by = "fgsea"
)


# --------------------------------------------------
# 6. Reactome GSEA
# --------------------------------------------------

set.seed(42)

cat("Running Reactome GSEA...\n")

gsea_reactome <- gsePathway(
    geneList = gene_list,
    organism = "human",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE
)


# --------------------------------------------------
# 7. Save + summarize
# --------------------------------------------------

save_gsea <- function(
    obj,
    prefix,
    title
) {

    df <- as.data.frame(obj)

    sig <- df[
        !is.na(df$p.adjust) &
        df$p.adjust < 0.05,
    ]

    write.table(
        df,
        paste0(
            "results/gsea/",
            prefix,
            "_all.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    write.table(
        sig,
        paste0(
            "results/gsea/",
            prefix,
            "_significant.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat("\n", prefix, "\n", sep="")
    cat("Significant pathways:",
        nrow(sig), "\n")

    cat("Positive NES (PD enriched):",
        sum(sig$NES > 0), "\n")

    cat("Negative NES (Healthy enriched):",
        sum(sig$NES < 0), "\n\n")

    if (nrow(sig) > 0) {

        cat("Top positive NES:\n")

        print(
            head(
                sig[
                    order(
                        sig$NES,
                        decreasing = TRUE
                    ),
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

        cat("\nTop negative NES:\n")

        print(
            head(
                sig[
                    order(
                        sig$NES
                    ),
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

        pdf(
            paste0(
                "results/gsea/",
                prefix,
                "_dotplot.pdf"
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
}


save_gsea(
    gsea_go,
    "GSEA_GO_BP",
    "GSEA - GO Biological Process"
)

save_gsea(
    gsea_reactome,
    "GSEA_Reactome",
    "GSEA - Reactome"
)


# --------------------------------------------------
# 8. Save objects
# --------------------------------------------------

saveRDS(
    gsea_go,
    "results/gsea/GSEA_GO_BP.rds"
)

saveRDS(
    gsea_reactome,
    "results/gsea/GSEA_Reactome.rds"
)

cat("\n=== GSEA completed ===\n")

