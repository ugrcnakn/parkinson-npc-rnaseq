
library(clusterProfiler)
library(ReactomePA)
library(enrichplot)
library(org.Hs.eg.db)

cat("=== Core pathway genes ===\n\n")

dir.create(
    "results/core_pathways",
    recursive = TRUE,
    showWarnings = FALSE
)

# --------------------------------------------------
# Load GSEA objects
# --------------------------------------------------

go <- readRDS(
    "results/gsea/GSEA_GO_BP.rds"
)

reactome <- readRDS(
    "results/gsea/GSEA_Reactome.rds"
)

ranked <- read.delim(
    "results/gsea/GSEA_ranked_genes.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# --------------------------------------------------
# Pathways of interest
# --------------------------------------------------

selected_reactome <- c(
    "R-HSA-6791226",
    "R-HSA-72613",
    "R-HSA-5620924",
    "R-HSA-5617833",
    "R-HSA-3238698"
)

selected_go <- c(
    "GO:0002181",
    "GO:0042254",
    "GO:0003341",
    "GO:0035082"
)

# --------------------------------------------------
# Extract leading-edge genes
# --------------------------------------------------

extract_core <- function(obj, selected_ids, database) {

    df <- as.data.frame(obj)

    df <- df[
        df$ID %in% selected_ids,
        ,
        drop = FALSE
    ]

    output <- list()

    for (i in seq_len(nrow(df))) {

        ids <- unlist(
            strsplit(
                df$core_enrichment[i],
                "/",
                fixed = TRUE
            )
        )

        ids <- unique(ids)

        map <- ranked[
            ranked$ENTREZID %in% ids,
            ,
            drop = FALSE
        ]

        map$Pathway_ID <- df$ID[i]
        map$Pathway <- df$Description[i]
        map$NES <- df$NES[i]
        map$p_adjust <- df$p.adjust[i]

        map$Direction <- ifelse(
            df$NES[i] > 0,
            "PD_enriched",
            "Healthy_enriched"
        )

        map$Database <- database

        output[[i]] <- map
    }

    do.call(rbind, output)
}

reactome_core <- extract_core(
    reactome,
    selected_reactome,
    "Reactome"
)

go_core <- extract_core(
    go,
    selected_go,
    "GO_BP"
)

core <- rbind(
    reactome_core,
    go_core
)

write.table(
    core,
    "results/core_pathways/core_pathway_genes.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# --------------------------------------------------
# Pathway summary
# --------------------------------------------------

summary_df <- unique(
    core[
        ,
        c(
            "Database",
            "Pathway_ID",
            "Pathway",
            "NES",
            "p_adjust",
            "Direction"
        )
    ]
)

write.table(
    summary_df,
    "results/core_pathways/core_pathway_summary.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("Selected pathways:\n\n")

print(summary_df)

# --------------------------------------------------
# Number of leading-edge genes
# --------------------------------------------------

cat("\nLeading-edge gene counts:\n\n")

counts <- table(
    core$Pathway
)

print(counts)

# --------------------------------------------------
# Print genes
# --------------------------------------------------

cat("\nCore genes by pathway:\n")

for (pathway in unique(core$Pathway)) {

    x <- core[
        core$Pathway == pathway,
        ,
        drop = FALSE
    ]

    cat(
        "\n--- ",
        pathway,
        " ---\n",
        sep = ""
    )

    genes <- unique(
        x$SYMBOL
    )

    genes <- genes[
        !is.na(genes) &
        genes != ""
    ]

    cat(
        paste(
            genes,
            collapse = ", "
        ),
        "\n"
    )
}

# --------------------------------------------------
# GSEA plots - Reactome
# --------------------------------------------------

pdf(
    "results/core_pathways/key_Reactome_GSEA_plots.pdf",
    width = 9,
    height = 6
)

reactome_df <- as.data.frame(reactome)

for (id in selected_reactome) {

    if (id %in% reactome_df$ID) {

        title <- reactome_df$Description[
            reactome_df$ID == id
        ][1]

        print(
            gseaplot2(
                reactome,
                geneSetID = id,
                title = title
            )
        )
    }
}

dev.off()

# --------------------------------------------------
# GSEA plots - GO
# --------------------------------------------------

pdf(
    "results/core_pathways/key_GO_GSEA_plots.pdf",
    width = 9,
    height = 6
)

go_df <- as.data.frame(go)

for (id in selected_go) {

    if (id %in% go_df$ID) {

        title <- go_df$Description[
            go_df$ID == id
        ][1]

        print(
            gseaplot2(
                go,
                geneSetID = id,
                title = title
            )
        )
    }
}

dev.off()

cat("\n=== Core pathway analysis completed ===\n")

