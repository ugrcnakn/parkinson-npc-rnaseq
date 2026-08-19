
library(DESeq2)

vsd <- readRDS("results/deseq2/vsd.rds")

mat <- assay(vsd)
meta <- as.data.frame(colData(vsd))

cor_matrix <- cor(mat, method = "pearson")
dist_matrix <- as.matrix(dist(t(mat)))

cat("=== 2B pairwise correlations ===\n\n")

samples_2B <- rownames(meta)[meta$line_code == "2B"]

print(
    round(
        cor_matrix[samples_2B, samples_2B],
        4
    )
)

cat("\n=== 2B pairwise distances ===\n\n")

print(
    round(
        dist_matrix[samples_2B, samples_2B],
        2
    )
)

cat("\n=== SRR16117933 correlations with all samples ===\n\n")

x <- sort(
    cor_matrix["SRR16117933", ],
    decreasing = TRUE
)

print(round(x, 4))

cat("\n=== SRR16117933 distances to all samples ===\n\n")

d <- sort(
    dist_matrix["SRR16117933", ]
)

print(round(d, 2))

