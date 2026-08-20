
OUT <- "results/enrichment"
SUMOUT <- "results/enrichment/summary"

dir.create(
    SUMOUT,
    recursive = TRUE,
    showWarnings = FALSE
)

read_safe <- function(file) {

    if (!file.exists(file)) {
        return(data.frame())
    }

    x <- read.delim(
        file,
        check.names = FALSE
    )

    return(x)
}

# ============================================================
# LOAD
# ============================================================

go_up <- read_safe(
    file.path(OUT, "ORA_GO_BP_PD_up.tsv")
)

go_down <- read_safe(
    file.path(OUT, "ORA_GO_BP_PD_down.tsv")
)

react_up <- read_safe(
    file.path(OUT, "ORA_Reactome_PD_up.tsv")
)

react_down <- read_safe(
    file.path(OUT, "ORA_Reactome_PD_down.tsv")
)

gsea_go <- read_safe(
    file.path(OUT, "GSEA_GO_BP.tsv")
)

gsea_react <- read_safe(
    file.path(OUT, "GSEA_Reactome.tsv")
)

# ============================================================
# ORA TOP TABLE
# ============================================================

top_ora <- function(x, n = 15) {

    if (nrow(x) == 0) {
        return(x)
    }

    x <- x[
        order(x$p.adjust),
    ]

    cols <- c(
        "Description",
        "GeneRatio",
        "Count",
        "p.adjust"
    )

    cols <- cols[
        cols %in% colnames(x)
    ]

    head(
        x[, cols, drop = FALSE],
        n
    )
}

# ============================================================
# GSEA TOP TABLE
# Positive NES = PD-up program
# Negative NES = PD-down program
# ============================================================

top_gsea <- function(x, direction, n = 15) {

    if (nrow(x) == 0) {
        return(x)
    }

    if (direction == "up") {

        x <- x[
            !is.na(x$NES) &
            x$NES > 0,
        ]

    } else {

        x <- x[
            !is.na(x$NES) &
            x$NES < 0,
        ]
    }

    x <- x[
        order(x$p.adjust),
    ]

    cols <- c(
        "Description",
        "NES",
        "setSize",
        "p.adjust"
    )

    cols <- cols[
        cols %in% colnames(x)
    ]

    head(
        x[, cols, drop = FALSE],
        n
    )
}

# ============================================================
# CREATE TABLES
# ============================================================

tables <- list(

    ORA_GO_UP =
        top_ora(go_up),

    ORA_GO_DOWN =
        top_ora(go_down),

    ORA_REACTOME_UP =
        top_ora(react_up),

    ORA_REACTOME_DOWN =
        top_ora(react_down),

    GSEA_GO_UP =
        top_gsea(gsea_go, "up"),

    GSEA_GO_DOWN =
        top_gsea(gsea_go, "down"),

    GSEA_REACTOME_UP =
        top_gsea(gsea_react, "up"),

    GSEA_REACTOME_DOWN =
        top_gsea(gsea_react, "down")
)

# ============================================================
# SAVE + PRINT
# ============================================================

for (name in names(tables)) {

    x <- tables[[name]]

    write.table(
        x,
        file.path(
            SUMOUT,
            paste0(name, "_top15.tsv")
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    cat("\n")
    cat("============================================\n")
    cat(name, "\n")
    cat("============================================\n")

    if (nrow(x) == 0) {

        cat("No significant terms.\n")

    } else {

        print(
            x,
            row.names = FALSE
        )
    }
}

