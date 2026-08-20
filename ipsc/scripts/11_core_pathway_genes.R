library(org.Hs.eg.db)
library(AnnotationDbi)

OUT <- "results/enrichment"
CORE <- "results/enrichment/core_genes"

dir.create(
    CORE,
    recursive = TRUE,
    showWarnings = FALSE
)

# ------------------------------------------------------------
# GSEA results
# ------------------------------------------------------------

go <- read.delim(
    file.path(OUT, "GSEA_GO_BP.tsv"),
    check.names = FALSE
)

react <- read.delim(
    file.path(OUT, "GSEA_Reactome.tsv"),
    check.names = FALSE
)

# ------------------------------------------------------------
# Themes we want to inspect
# ------------------------------------------------------------

themes <- list(

    translation = c(
        "ribosome",
        "translation",
        "rRNA",
        "ribonucleoprotein"
    ),

    cilia = c(
        "cilium",
        "ciliary",
        "basal body",
        "intraciliary"
    ),

    microtubule_centrosome = c(
        "microtubule",
        "centrosome",
        "spindle"
    ),

    development = c(
        "morphogenesis",
        "development",
        "differentiation"
    )
)

# ------------------------------------------------------------
# Parse core_enrichment Entrez IDs
# ------------------------------------------------------------

extract_theme <- function(df, patterns, source) {

    pattern <- paste(patterns, collapse = "|")

    sub <- df[
        grepl(
            pattern,
            df$Description,
            ignore.case = TRUE
        ),
    ]

    if (nrow(sub) == 0) {
        return(data.frame())
    }

    rows <- list()

    for (i in seq_len(nrow(sub))) {

        ids <- unlist(
            strsplit(
                as.character(sub$core_enrichment[i]),
                "/"
            )
        )

        ids <- ids[ids != ""]

        if (length(ids) == 0) {
            next
        }

        symbols <- mapIds(
            org.Hs.eg.db,
            keys = ids,
            column = "SYMBOL",
            keytype = "ENTREZID",
            multiVals = "first"
        )

        tmp <- data.frame(
            source = source,
            pathway = sub$Description[i],
            NES = sub$NES[i],
            padj = sub$p.adjust[i],
            entrez = ids,
            symbol = as.character(symbols),
            stringsAsFactors = FALSE
        )

        rows[[length(rows) + 1]] <- tmp
    }

    if (length(rows) == 0) {
        return(data.frame())
    }

    do.call(
        rbind,
        rows
    )
}

# ------------------------------------------------------------
# Export each theme
# ------------------------------------------------------------

for (theme in names(themes)) {

    go_theme <- extract_theme(
        go,
        themes[[theme]],
        "GO_BP"
    )

    react_theme <- extract_theme(
        react,
        themes[[theme]],
        "Reactome"
    )

    combined <- rbind(
        go_theme,
        react_theme
    )

    write.table(
        combined,
        file.path(
            CORE,
            paste0(theme, "_leading_edge.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat("\n")
    cat("====================================\n")
    cat(toupper(theme), "\n")
    cat("====================================\n")

    if (nrow(combined) == 0) {
        cat("No genes found.\n")
        next
    }

    counts <- sort(
        table(combined$symbol),
        decreasing = TRUE
    )

    counts <- counts[
        names(counts) != "NA"
    ]

    print(
        head(counts, 25)
    )
}

