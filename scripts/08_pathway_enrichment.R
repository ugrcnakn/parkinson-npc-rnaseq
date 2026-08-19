
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

cat("=== Robust DEG pathway enrichment ===\n\n")

dir.create(
    "results/enrichment",
    recursive = TRUE,
    showWarnings = FALSE
)

# --------------------------------------------------
# 1. Load data
# --------------------------------------------------

robust <- read.delim(
    "results/deseq2/robust_1099_DEGs_annotated.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

all_tested <- read.delim(
    "results/deseq2/all12_DESeq2_results_annotated.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# --------------------------------------------------
# 2. Protein-coding genes only
# --------------------------------------------------

robust_pc <- robust[
    robust$gene_biotype == "protein_coding",
]

background_pc <- all_tested[
    all_tested$gene_biotype == "protein_coding",
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

cat("Robust protein-coding DEG:",
    nrow(robust_pc), "\n")

cat("PD-up Ensembl genes:",
    length(up_ens), "\n")

cat("PD-down Ensembl genes:",
    length(down_ens), "\n")

cat("Background protein-coding genes:",
    length(background_ens), "\n\n")


# --------------------------------------------------
# 3. ENSEMBL -> ENTREZ conversion
# --------------------------------------------------

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
    "results/enrichment/background_id_mapping.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    up_map,
    "results/enrichment/PD_up_id_mapping.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    down_map,
    "results/enrichment/PD_down_id_mapping.tsv",
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

cat("Mapped background Entrez:",
    length(background_entrez), "\n")

cat("Mapped PD-up Entrez:",
    length(up_entrez), "\n")

cat("Mapped PD-down Entrez:",
    length(down_entrez), "\n\n")


# --------------------------------------------------
# 4. GO Biological Process enrichment
# --------------------------------------------------

cat("Running GO BP - PD up...\n")

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

cat("Running GO BP - PD down...\n")

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


# --------------------------------------------------
# 5. Reactome enrichment
# --------------------------------------------------

cat("Running Reactome - PD up...\n")

reactome_up <- enrichPathway(
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

cat("Running Reactome - PD down...\n")

reactome_down <- enrichPathway(
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


# --------------------------------------------------
# 6. Save function
# --------------------------------------------------

save_result <- function(
    obj,
    prefix,
    title
) {

    df <- as.data.frame(obj)

    write.table(
        df,
        paste0(
            "results/enrichment/",
            prefix,
            ".tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat(
        prefix,
        "- significant terms:",
        nrow(df),
        "\n"
    )

    if (nrow(df) > 0) {

        cat("\nTop terms for", prefix, ":\n")

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
            paste0(
                "results/enrichment/",
                prefix,
                "_dotplot.pdf"
            ),
            width = 10,
            height = 8
        )

        print(
            enrichplot::dotplot(
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


# --------------------------------------------------
# 7. Save all results
# --------------------------------------------------

save_result(
    go_up,
    "GO_BP_PD_up",
    "GO Biological Process - Up in PD"
)

save_result(
    go_down,
    "GO_BP_PD_down",
    "GO Biological Process - Down in PD"
)

save_result(
    reactome_up,
    "Reactome_PD_up",
    "Reactome pathways - Up in PD"
)

save_result(
    reactome_down,
    "Reactome_PD_down",
    "Reactome pathways - Down in PD"
)


# --------------------------------------------------
# 8. Save R objects
# --------------------------------------------------

saveRDS(
    go_up,
    "results/enrichment/GO_BP_PD_up.rds"
)

saveRDS(
    go_down,
    "results/enrichment/GO_BP_PD_down.rds"
)

saveRDS(
    reactome_up,
    "results/enrichment/Reactome_PD_up.rds"
)

saveRDS(
    reactome_down,
    "results/enrichment/Reactome_PD_down.rds"
)

cat("=== Enrichment analysis completed ===\n")

