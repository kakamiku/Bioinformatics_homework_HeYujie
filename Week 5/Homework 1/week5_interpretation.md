# Week 5 Homework 1 — Interpretation

**Comparison:** treated versus control, design `~ batch + condition`, reference level `control`,
coefficient `condition_treated_vs_control`.

**Thresholds reported:** statistical threshold `padj < 0.05` (Benjamini–Hochberg FDR) **and**
effect-size threshold `|log2FoldChange| >= 1`. The effect-size cut is applied to the
**apeglm-shrunken** estimate. For transparency: the FDR cut alone gives 83 genes, FDR plus the
*unshrunken* `|log2FC| >= 1` gives 74, and the reported 60 is the shrunken combination.

---

## Interpretation (147 words)

Treated versus control was compared under ~ batch + condition (reference: control); 989 of
1,000 genes passed the 10-counts-in-3-samples filter. The strongest QC observation is the VST
PCA: PC1 (24% variance) separates all six controls from all six treated samples with no overlap
and shows no batch signal, while batch gives only a weak, incomplete grouping on PC2 (9%), so
batch belongs in the model rather than in a sample-removal decision. At padj below 0.05 and
shrunken absolute log2FC of at least 1, 60 genes are significant: 36 up, 24 down in treated. Only
~6% of tested genes respond, so the treatment acts on a restricted target set, not a global
reprogramming. Limitation: two replicates per batch-condition cell, and synthetic gene
identifiers, so no pathway can be named. AI drafted the direction-check block; I re-derived the
fold change by hand and verified the coefficient name against resultsNames.

---

## Supporting numbers (reference table, not part of the 147-word paragraph)

| Quantity | Value |
|---|---|
| Genes before / after filtering | 1,000 / 989 |
| `resultsNames(dds)` | Intercept, batch_B_vs_A, batch_C_vs_A, condition_treated_vs_control |
| PC1 / PC2 variance explained | 24% / 9% |
| PC1 by condition (one-way ANOVA) | p = 4.1e-12 — complete separation |
| PC1 by batch | p = 0.99 — no batch signal on PC1 |
| PC2 by batch | p = 0.048 — weak; PC2 is dominated by within-batch scatter (TB1 and TB2 are the two most extreme PC2 values and are the same batch and condition) |
| Significant (padj < 0.05 and abs shrunken log2FC >= 1) | 60 — 36 up, 24 down |
| Top gene by effect size | Gene0008, shrunken log2FC 1.93 (up in treated) |
| Library sizes | 112,619 – 151,748 (≈35% spread) |
