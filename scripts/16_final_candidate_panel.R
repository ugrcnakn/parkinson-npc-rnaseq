
cat("=== FINAL CANDIDATE GENE PANEL ===\n\n")

res <- read.delim(
    "results/deseq2_paired/paired_DESeq2_results_annotated.tsv",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

genes <- c(

    # WNT
    "AXIN2",
    "RSPO1",
    "APCDD1",
    "DKK1",
    "WNT1",
    "WNT2B",
    "WNT3A",
    "WNT4",
    "WNT9A",
    "WNT10B",

    # Cilia / IFT
    "FOXJ1",
    "RFX3",
    "IFT88",
    "IFT172",
    "IFT81",
    "IFT56",
    "IFT46",
    "IFT74",
    "KIF3B",
    "BBS4",
    "CEP290",

    # Neuronal / NMDA
    "GRIN2B",
    "GRIN2C",
    "GRIN3A",
    "DLG2",
    "ERBB4",

    # Translation / ribosome
    "RPL22L1",
    "DDX21",
    "NPM1",
    "EIF4E",
    "EIF4G1",
    "RPS6",

    # Metabolic / top DEG
    "BCAT1",
    "HK2"
)

panel <- res[
    res$gene_name %in% genes,
    c(
        "Geneid",
        "gene_name",
        "gene_biotype",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj"
    )
]

panel <- panel[
    match(
        genes,
        panel$gene_name
    ),
]

panel <- panel[
    !is.na(panel$gene_name),
]

panel$direction <- ifelse(
    panel$log2FoldChange > 0,
    "Higher_in_PD",
    "Lower_in_PD"
)

write.table(
    panel,
    "results/final_core_biology/final_candidate_gene_panel.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

print(
    panel,
    row.names = FALSE
)

cat(
    "\n=== FINAL CANDIDATE PANEL COMPLETED ===\n"
)

