# Parkinson's Disease Neural Progenitor Cell RNA-seq Analysis

## Overview

This project performs an end-to-end bulk RNA-seq analysis of human neural
progenitor cells (NPCs) derived from monozygotic twin pairs discordant for
Parkinson's disease (PD).

The analysis starts from raw FASTQ files and proceeds through read alignment,
gene-level quantification, differential expression analysis, pathway
enrichment, sensitivity analysis, and biological interpretation.

The final statistical model accounts for the paired twin structure:

    design = ~ pair + condition

with Healthy used as the reference condition.

---

## Experimental Design

The dataset contains 12 single-end RNA-seq samples.

Two monozygotic twin pairs were analyzed:

- Pair 1
  - Healthy (1K): 3 RNA-seq replicates
  - PD (1B): 3 RNA-seq replicates

- Pair 2
  - Healthy (2K): 3 RNA-seq replicates
  - PD (2B): 3 RNA-seq replicates

Total:

- 6 PD samples
- 6 Healthy samples

---

## Workflow

Raw FASTQ

↓

Read alignment with HISAT2

↓

Sorted and indexed BAM files

↓

Gene-level quantification with featureCounts

↓

Count matrix construction

↓

DESeq2 normalization and QC

↓

Paired differential expression analysis

↓

Outlier sensitivity analysis

↓

Gene annotation

↓

Volcano plot and DEG heatmap

↓

GO Biological Process enrichment

↓

Reactome pathway enrichment

↓

Gene Set Enrichment Analysis (GSEA)

↓

Leading-edge gene analysis

↓

Biological interpretation

---

## Reference Genome

Genome:

- Human GRCh38

Annotation:

- Ensembl release 92

Alignment:

- HISAT2

Quantification:

- featureCounts

Differential expression:

- DESeq2

Functional enrichment:

- clusterProfiler
- ReactomePA

---

## Alignment

Most samples showed high alignment rates.

Observed alignment rates ranged approximately from 80% to 98%.

One sample, SRR16117933, showed unusual behavior in PCA and sample-distance
analysis and was therefore investigated with a sensitivity analysis rather
than automatically removed.

---

## Differential Expression

Genes were retained when:

    count >= 10 in at least 6 samples

Genes tested:

    15,373

Differential expression threshold:

    adjusted p-value < 0.05
    |log2 fold change| >= 1

Final paired DESeq2 analysis:

    Significant DEGs: 1,207

    Higher in PD: 682
    Lower in PD: 525

---

## Sensitivity Analysis

The paired analysis was repeated after removing SRR16117933.

Results:

    Paired all-sample DEGs: 1,207
    Without SRR16117933: 1,335
    Common significant DEGs: 1,066

    Global log2FC correlation: 0.955
    Direction concordance among common DEGs: 100%

Because the global differential-expression pattern remained highly
consistent, SRR16117933 was retained in the primary analysis.

The 1,066 DEGs significant in both analyses were used as a robust DEG set.

---

## Robust Differentially Expressed Genes

Robust DEGs:

    1,066

Robust PD-up:

    614

Robust PD-down:

    452

Robust protein-coding DEGs:

    871

Protein-coding PD-up:

    532

Protein-coding PD-down:

    339

---

## Major Biological Findings

### 1. Reduced cilia-related transcriptional programs in PD NPCs

Both over-representation analysis and GSEA identified reduced activity of
cilia-associated pathways.

Major signals included:

- cilium movement
- axoneme assembly
- cilium assembly
- intraflagellar transport
- microtubule-associated movement

Representative genes included:

- FOXJ1
- BBS4
- IFT46
- IFT74
- IFT88
- IFT81
- IFT172
- KIF3B
- CEP290

These findings suggest coordinated alteration of cilia-associated
transcriptional programs in PD-derived NPCs.

---

### 2. Altered WNT-associated transcription

Several WNT-related genes showed lower expression in PD NPCs.

Examples:

- AXIN2
- RSPO1
- APCDD1
- DKK1
- WNT1
- WNT2B
- WNT3A
- WNT4
- WNT9A
- WNT10B

Reactome enrichment also identified:

- WNT ligand biogenesis and trafficking
- regulation of TCF-dependent WNT signaling

as enriched among genes lower in PD.

---

### 3. Reduced neuronal signaling programs

Pathway analysis identified decreased enrichment of neuronal processes,
including:

- Neuronal System
- NMDA receptor assembly and cell-surface presentation

Representative genes included:

- GRIN2B
- GRIN2C
- GRIN3A
- DLG2
- ERBB4

---

### 4. Increased ribosome and translation-associated programs

GSEA identified strong positive enrichment in PD for:

- ribosome biogenesis
- rRNA processing
- cytoplasmic translation
- translation initiation
- ribonucleoprotein complex biogenesis

Representative genes included:

- RPL22L1
- DDX21
- NPM1
- EIF4E
- EIF4G1

This signal should be interpreted as a transcriptomic program specific to
this NPC model rather than a universal increase in translation in
Parkinson's disease.

---

## Selected Candidate Genes

Examples of strongly altered genes include:

| Gene | log2FC | Direction |
|------|-------:|-----------|
| AXIN2 | -2.29 | Lower in PD |
| RSPO1 | -3.57 | Lower in PD |
| DKK1 | -3.07 | Lower in PD |
| WNT3A | -2.14 | Lower in PD |
| WNT4 | -3.80 | Lower in PD |
| FOXJ1 | -2.39 | Lower in PD |
| GRIN2B | -1.40 | Lower in PD |
| ERBB4 | -3.62 | Lower in PD |
| RPL22L1 | +1.27 | Higher in PD |
| DDX21 | +1.05 | Higher in PD |
| BCAT1 | +1.43 | Higher in PD |
| HK2 | +1.44 | Higher in PD |

Positive log2FC indicates higher expression in PD.

Negative log2FC indicates lower expression in PD.

---

## DKK1

DKK1 remained significant in the paired analysis:

    log2FC = -3.07
    adjusted p-value = 0.00137

The direction of the DKK1 change was also preserved in the sensitivity
analysis excluding SRR16117933.

---

## Low-Expression Candidate Genes

Some genes highlighted in the original study were not included in the final
DESeq2 analysis because they did not pass the predefined expression filter.

Examples:

WNT7A:

    total raw counts = 32
    count >= 10 in 1/12 samples

TNF:

    total raw counts = 87
    count >= 10 in 3/12 samples

INHBA:

    total raw counts = 233
    count >= 10 in 4/12 samples

These genes were not forced into the differential-expression analysis.

---

## Interpretation

The combined differential-expression and pathway analyses suggest a
transcriptional phenotype in PD-derived neural progenitor cells characterized
by:

1. reduced cilia and intraflagellar transport programs,
2. altered WNT-associated developmental signaling,
3. reduced neuronal signaling programs,
4. increased ribosome, rRNA-processing, and translation-associated programs.

These observations represent associations within this experimental model and
should not be interpreted as proof of causal mechanisms in Parkinson's
disease.

---

## Limitations

This dataset contains two monozygotic twin pairs with three RNA-seq
replicates per condition within each pair.

Although the paired statistical design controls for the twin-pair effect,
the number of independent biological pairs is small.

Technical or culture replicates should not be interpreted as independent
patients.

Therefore, the results should be considered exploratory and require
validation in larger independent cohorts and functional experiments.

The analysis is based on bulk RNA-seq and therefore cannot distinguish
cell-to-cell heterogeneity.

---

## Key Output Directories

    results/alignment/
    results/counts/
    results/deseq2_paired/
    results/deseq2_paired/final/
    results/enrichment_paired_final/
    results/final_core_biology/

---

## Key Final Files

    paired_DESeq2_results_annotated.tsv
    paired_robust_1066_DEGs.tsv
    paired_volcano_PD_vs_Healthy.pdf
    paired_top30_robust_DEG_heatmap.pdf

    ORA_GO_BP_PD_up.tsv
    ORA_GO_BP_PD_down.tsv

    ORA_Reactome_PD_up.tsv
    ORA_Reactome_PD_down.tsv

    GSEA_GO_BP_paired_significant.tsv
    GSEA_Reactome_paired_significant.tsv

    final_leading_edge_genes.tsv
    final_candidate_gene_panel.tsv

---

## Project Status

End-to-end RNA-seq analysis completed.

Raw FASTQ -> alignment -> quantification -> paired differential expression ->
sensitivity analysis -> enrichment -> GSEA -> biological interpretation.

