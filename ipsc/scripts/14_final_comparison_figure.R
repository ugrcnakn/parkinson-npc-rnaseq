library(ggplot2)

OUT <- "results/final_figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. ALL SHARED TESTED GENES
# ============================================================

all_genes <- read.delim(
    "results/ipsc_vs_npc/all_shared_tested_genes_iPSC_NPC.tsv",
    check.names = FALSE
)

all_genes$group <- "Other"

all_genes$group[
    abs(all_genes$iPSC_log2FC) >= 0.5 &
    abs(all_genes$NPC_log2FC) >= 0.5
] <- "LFC>=0.5 in both"

p1 <- ggplot(
    all_genes,
    aes(
        x = iPSC_log2FC,
        y = NPC_log2FC
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dotted"
    ) +
    geom_point(
        aes(alpha = group),
        size = 1
    ) +
    scale_alpha_manual(
        values = c(
            "Other" = 0.15,
            "LFC>=0.5 in both" = 0.75
        )
    ) +
    theme_bw() +
    guides(alpha = "none") +
    labs(
        title = "Global PD effect: iPSC vs NPC",
        subtitle = "Shared genes tested in both cell states",
        x = "iPSC log2 fold change",
        y = "NPC log2 fold change"
    )

ggsave(
    file.path(
        OUT,
        "global_iPSC_vs_NPC_log2FC.png"
    ),
    p1,
    width = 8,
    height = 7,
    dpi = 300
)

ggsave(
    file.path(
        OUT,
        "global_iPSC_vs_NPC_log2FC.pdf"
    ),
    p1,
    width = 8,
    height = 7
)

# ============================================================
# 2. SHARED ROBUST 50
# ============================================================

robust <- read.delim(
    "results/ipsc_vs_npc/shared_50_robust_genes_iPSC_NPC.tsv",
    check.names = FALSE
)

p2 <- ggplot(
    robust,
    aes(
        x = iPSC_log2FC,
        y = NPC_log2FC
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dotted"
    ) +
    geom_point(
        aes(shape = direction),
        size = 3,
        alpha = 0.8
    ) +
    theme_bw() +
    labs(
        title = "Shared robust PD-associated genes",
        subtitle = "50 genes common to iPSC and NPC robust sets",
        x = "iPSC log2 fold change",
        y = "NPC log2 fold change",
        shape = "Direction"
    )

ggsave(
    file.path(
        OUT,
        "shared_50_robust_iPSC_vs_NPC.png"
    ),
    p2,
    width = 8,
    height = 7,
    dpi = 300
)

ggsave(
    file.path(
        OUT,
        "shared_50_robust_iPSC_vs_NPC.pdf"
    ),
    p2,
    width = 8,
    height = 7
)

cat("\n")
cat("========================================\n")
cat("FINAL FIGURES CREATED\n")
cat("========================================\n")
cat("global_iPSC_vs_NPC_log2FC.png\n")
cat("global_iPSC_vs_NPC_log2FC.pdf\n")
cat("shared_50_robust_iPSC_vs_NPC.png\n")
cat("shared_50_robust_iPSC_vs_NPC.pdf\n")
cat("========================================\n")
