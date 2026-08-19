
library(clusterProfiler)
library(org.Hs.eg.db)

cat("=== FINAL CORE BIOLOGY ===\n\n")

outdir <- "results/final_core_biology"

dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
)

# ==================================================
# 1. LOAD FINAL OBJECTS
# ==================================================

go <- readRDS(
    "results/enrichment_paired_final/GSEA_GO_BP_paired.rds"
)

react <- readRDS(
    "results/enrichment_paired_final/GSEA_Reactome_paired.rds"
)

ranked <- read.delim(
    "results/enrichment_paired_final/paired_GSEA_ranked_genes.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

full <- read.delim(
    "results/deseq2_paired/paired_DESeq2_results_annotated.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

# ==================================================
# 2. FINAL KEY GSEA PATHWAYS
# ==================================================

go_ids <- c(
    "GO:0042254",   # ribosome biogenesis
    "GO:0002181",   # cytoplasmic translation
    "GO:0035082",   # axoneme assembly
    "GO:0003341"    # cilium movement
)

react_ids <- c(
    "R-HSA-6791226", # rRNA processing
    "R-HSA-72613",   # translation initiation
    "R-HSA-5620924", # intraflagellar transport
    "R-HSA-5617833"  # cilium assembly
)

extract_core <- function(obj, ids, database) {

    df <- as.data.frame(obj)

    df <- df[
        df$ID %in% ids,
        ,
        drop = FALSE
    ]

    output <- list()

    for (i in seq_len(nrow(df))) {

        entrez <- unique(
            strsplit(
                df$core_enrichment[i],
                "/",
                fixed = TRUE
            )[[1]]
        )

        x <- ranked[
            as.character(ranked$ENTREZID) %in% entrez,
            ,
            drop = FALSE
        ]

        x$Database <- database
        x$Pathway_ID <- df$ID[i]
        x$Pathway <- df$Description[i]
        x$NES <- df$NES[i]
        x$p_adjust <- df$p.adjust[i]

        x$Direction <- ifelse(
            df$NES[i] > 0,
            "PD_enriched",
            "Healthy_enriched"
        )

        output[[i]] <- x
    }

    do.call(
        rbind,
        output
    )
}

go_core <- extract_core(
    go,
    go_ids,
    "GO_BP"
)

react_core <- extract_core(
    react,
    react_ids,
    "Reactome"
)

core <- rbind(
    go_core,
    react_core
)

write.table(
    core,
    file.path(
        outdir,
        "final_leading_edge_genes.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ==================================================
# 3. PATHWAY SUMMARY
# ==================================================

summary <- unique(
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

cat("FINAL PATHWAYS:\n\n")

print(summary)

cat("\nLeading-edge counts:\n\n")

print(
    table(core$Pathway)
)

# ==================================================
# 4. PRINT CORE GENES
# ==================================================

cat("\nCORE GENES:\n")

for (p in unique(core$Pathway)) {

    x <- core[
        core$Pathway == p,
        ,
        drop = FALSE
    ]

    genes <- unique(x$SYMBOL)

    genes <- genes[
        !is.na(genes) &
        genes != ""
    ]

    cat(
        "\n--- ",
        p,
        " ---\n",
        sep = ""
    )

    cat(
        paste(
            genes,
            collapse = ", "
        ),
        "\n"
    )
}

# ==================================================
# 5. WNT ORA RESULTS
# ==================================================

react_down <- read.delim(
    "results/enrichment_paired_final/ORA_Reactome_PD_down.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

wnt <- react_down[
    grepl(
        "WNT|Neuronal|NMDA",
        react_down$Description,
        ignore.case = TRUE
    ),
]

write.table(
    wnt,
    file.path(
        outdir,
        "WNT_neuronal_ORA.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat("\nWNT / neuronal ORA:\n\n")

print(
    wnt[
        ,
        c(
            "ID",
            "Description",
            "GeneRatio",
            "p.adjust",
            "geneID"
        )
    ]
)

# ==================================================
# 6. ORIGINAL PAPER CANDIDATES
# ==================================================

targets <- c(
    "TNF",
    "INHBA",
    "WNT7A",
    "DKK1"
)

cat("\nOriginal paper candidate genes in paired DESeq2:\n\n")

for (g in targets) {

    x <- full[
        full$gene_name == g,
        ,
        drop = FALSE
    ]

    cat("\n--- ", g, " ---\n", sep="")

    if (nrow(x) == 0) {

        cat("Not present after expression filtering.\n")

    } else {

        print(
            x[
                ,
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
    }
}

# ==================================================
# 7. RAW COUNT FILTER CHECK
# ==================================================

counts <- read.delim(
    "results/counts/count_matrix.tsv",
    row.names = 1,
    check.names = FALSE
)

anno <- read.delim(
    "results/annotation/ensembl92_gene_annotation.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

target_anno <- anno[
    anno$gene_name %in% targets,
]

cat("\nRaw-count filtering check:\n\n")

for (i in seq_len(nrow(target_anno))) {

    id <- target_anno$Geneid[i]
    gene <- target_anno$gene_name[i]

    cat(
        "\n",
        gene,
        " (",
        id,
        ")\n",
        sep = ""
    )

    if (!(id %in% rownames(counts))) {

        cat("Not present in count matrix.\n")
        next
    }

    values <- as.numeric(
        counts[id, ]
    )

    cat(
        "Total raw counts:",
        sum(values),
        "\n"
    )

    cat(
        "Samples with count >=10:",
        sum(values >= 10),
        "/12\n"
    )

    cat(
        "Passes our filter:",
        sum(values >= 10) >= 6,
        "\n"
    )

    cat(
        "Counts:",
        paste(
            values,
            collapse = ", "
        ),
        "\n"
    )
}

cat(
    "\n=== FINAL CORE BIOLOGY COMPLETED ===\n"
)

