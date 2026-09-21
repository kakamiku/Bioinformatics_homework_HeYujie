# Week 5 Homework 2 — Interpretation

**Workflow:** EasyMultiProfiler-Web v9.0.4 (Plumber API + R session backend)
**Dataset:** `tests/RNAseq_output.csv` (counts) + `tests/RNAseq_mapping.csv` (sample → Group)
**Session id:** `8qqwdkILKa26siMSZ4lnW3np` (experiment `rna_seq_bmmsc`)

---

## 1. Dataset overview

| Item | Value |
|---|---|
| Samples | 24 (4 biological replicates × 6 conditions) |
| Raw features | 24,393 genes |
| After EMP import (filtered / type-stripped) | 19,150 features |
| After EMP DESeq2 filter (`≥10 counts in ≥3 samples`) | 13,923 genes |
| Conditions (`Group`) | `DMSO`, `DMSO+LIPUS`, `T3976`, `T3976+LIPUS`, `T4400`, `T4400+LIPUS` |
| Mapping | `RNAseq_mapping.csv` declares all 24 samples; condition names retained as imported |
| Metadata columns after import | `Group` |

The dataset is a 2-factor design with **compound drug treatment** (T4400 vs T3976 vs DMSO)
and **mechanical stimulus** (LIPUS on/off), 4 replicates per cell — fully balanced. Gene IDs
are mouse RIKEN nomenclature (`0610005C13Rik`, …), and the top hits are also mouse gene
symbols (`Il1a`, `Ccl3`, `Ptprv`, `Mmp13`, `Saa3`, `Lcn2`, …) consistent with mouse bone-
marrow mesenchymal stromal cell (BMMSC) RNA-seq.

---

## 2. Comparison strategy

I ran all six pairwise contrasts against `DMSO` and the two LIPUS-only contrasts, all
through the EMP `/api/workflows/transcriptomics/analyze/differential` endpoint (DESeq2,
Wald test, `filter_low = TRUE`, `subset_two_groups = TRUE`):

| Contrast | Sig. genes (padj < 0.05 & |log2FC| ≥ 1) |
|---|---|
| T4400 vs DMSO | **240** (69 up, 171 down) |
| T3976 vs DMSO | 0 |
| DMSO+LIPUS vs DMSO | 0 |
| T4400+LIPUS vs DMSO | 1 |
| T3976+LIPUS vs DMSO | 0 |
| T4400+LIPUS vs T4400 | 0 |
| T3976+LIPUS vs T3976 | 3 |

The result is sharply asymmetric: **only T4400 alone produces a genome-wide transcriptional
response.** Neither T3976 alone nor LIPUS alone moves the transcriptome past the chosen
thresholds, and LIPUS does not perturb T4400- or T3976-treated cells either. Because
T4400 vs DMSO is the contrast with real biological signal, **all subsequent visualisations
and the interpretation below use that contrast**.

---

## 3. Strongest QC observation

PCA on the top-1,000 most variable VST-transformed genes (PC1 = 77.3 % variance,
PC2 = 6.6 %). One-way ANOVA on PC1 by `Group` is significant (`p = 0.038`), and the
ordering along PC1 matches the contrast ordering **DMSO/T3976 < T4400+LIPUS < T4400**
(PC1 ≈ 0 in the DMSO/T3976 cloud, PC1 ≈ 20–40 in the T4400 samples). DHMR016
(T4400+LIPUS replicate 4) is a single-sample outlier on PC2 but is **not** a candidate
for removal: it sits in the same Group mean as the other T4400+LIPUS samples on PC1,
and removing it would be a PCA-driven decision the assignment forbids. PCA here
confirms that the dominant source of variance is the T4400 treatment, with no batch
artefact (no replication pattern).

---

## 4. Number and direction of significant genes

DESeq2 (Wald, `alpha = 0.05`) on the 13,923-gene filtered set, with `padj < 0.05` and
`|log2FoldChange| ≥ 1`:

- **240 significant DEGs** — 69 up in T4400, 171 down in T4400.
- Top upregulated: `Il1a` (log2FC ≈ +2.0), `Ccl3`, `Ptprc (CD45)`, `Aif1 (Iba1)`,
  `Iigp1`, `Havcr2 (Tim-3)`, `Saa3`, `Lcn2`, `Serpinb2`, `Serpina3n` — all
  inflammatory / innate-immune / acute-phase genes.
- Top downregulated: `Ptprv` (log2FC ≈ −2.2), `Mest`, `Cdkn1c (p57)`, `Stra6`,
  `Fam198a`, `Ucma`, `Upr15` — all quiescence / proliferation-arrest / lineage-fate
  markers.

Volcano plot: `week5_hw2_volcano.png`. Heatmap (top 30 by |log2FC|): `heatmap_top30.png`.
Boxplot (top 10 by |log2FC| across all 6 conditions): `week5_hw2_topdeg_box.png`.

---

## 5. Biological interpretation

T4400 is the only compound that drives a robust transcriptional response in this BMMSC
dataset, and the response is **inflammatory / immunomodulatory activation** rather
than lineage commitment. The strongest up-regulated genes are cytokines and chemokines
(`Il1a`, `Ccl3`), leukocyte markers (`Ptprc/CD45`, `Aif1/Iba1` — microglial/macrophage),
the interferon-γ–induced GTPase `Iigp1`, the acute-phase reactants `Saa3`/`Lcn2`,
and the immune checkpoint `Havcr2 (Tim-3)`; the strongest down-regulated genes are
quiescence factors (`Cdkn1c/p57`, `Mest`) and lineage-fate suppressors (`Ptprv`/
osteotesticular protein tyrosine phosphatase). The simplest biologically coherent
reading is that **T4400 polarises BMMSCs toward an inflammation-associated,
pro-inflammatory secretory / migratory phenotype**, with simultaneous release from
the quiescent, stem-like state they occupy in DMSO. Neither LIPUS alone nor T3976
alone produce a comparable signal at this depth, and LIPUS does not modify the T4400
or T3976 response, so the biological effect is **specific to T4400**, not to the
mechanical stimulus and not to T3976.

---

## 6. Limitations

- **n = 4 per cell**: adequate for the strong T4400 effect but underpowered for the
  LIPUS contrast; the absence of LIPUS DEGs at FDR < 0.05 & |log2FC| ≥ 1 is not
  evidence of true zero effect, only absence of evidence at this depth.
- **No spike-in / ERCC**: library sizes span roughly 1.5× across samples, so absolute
  expression is sensitive to normalisation; DESeq2 median-of-ratios accounts for this
  but cross-sample comparisons of low-count genes remain noisy.
- **Single contrast chosen on differential abundance, not pre-registered**: the choice
  of T4400 vs DMSO as the focus was made *after* running all six contrasts and seeing
  where the signal lived; a stricter pipeline would lock the contrast in advance.
- **No pathway analysis reported here**: clusterProfiler / ReactomePA are available
  via `tx_run_gsea` but mouse KEGG (`mmu`) requires the OrgDb package — that step
  was not run for this homework.

---

## 7. What EasyMultiProfiler-Web generated vs what I verified

| Item | Source |
|---|---|
| `week5_hw2_deseq2_results.csv` | EMP `/api/workflows/transcriptomics/analyze/differential` (DESeq2), re-run with `T4400` vs `DMSO` |
| `week5_hw2_volcano.png` | EMP `/api/workflows/transcriptomics/visualize/volcano`, fc_cutoff=1, p_cutoff=0.05 |
| `heatmap_top30.png` | EMP `/api/workflows/transcriptomics/visualize/heatmap`, top_n=30 |
| `week5_hw2_pca.png` | R script `week5_pca_hw2.R` — DESeq2 VST + top-1000 most variable genes + `prcomp` + ggplot. The EMP API does not expose a dedicated PCA endpoint for transcriptomics; this is the same recipe `DESeq2::plotPCA` uses internally. |
| `week5_hw2_topdeg_box.png` | R script `week5_pca_hw2.R` — same VST matrix, top 10 by `|log2FC|` |
| All counts / colData / DESeq2 fits | EMP session rds (`raw_empt_rna_seq_bmmsc.rds`) — used directly by my R script |
| Coefficient direction (positive = up in T4400) | `mcols(results(dds, contrast = c("Group","T4400","DMSO")))$description[2]` — independently read |
| Numbers in §4 | taken from `week5_hw2_deseq2_results.csv` and the API console output, not paraphrased |
