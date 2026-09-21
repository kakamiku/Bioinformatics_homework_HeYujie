# Week 5 Homework 1 — AI Use and Verification Log

**Assistant used:** Claude (Anthropic), `claude-opus-5`
**Date:** 2026-09-21
**Scope of AI use:** one task only — drafting the *independent direction-verification block*
(Section 7b of `week5_deseq2_analysis.R`). Everything else — the design choice, the filtering
rule, the thresholds, the figures, and the interpretation — was written and decided by me,
starting from the provided `Week5_Homework_Starter.R`. No other section of the script was
AI-generated.

---

## 1. The task given to AI

I wanted a second, model-independent way to confirm that the reported contrast really is
*treated relative to control* and not the reverse — i.e. that a positive `log2FoldChange` means
"higher in treated". DESeq2 will happily return a sign-flipped result if `relevel()` or the
`contrast` argument is wrong, and the sign is exactly the thing a reader will trust without
checking.

## 2. Prompt, preserved verbatim

> I have a DESeq2 object `dds` fitted with `design = ~ batch + condition`, where `condition` has
> levels `control` (reference) and `treated`. I extracted results with
> `results(dds, contrast = c("condition","treated","control"))` and shrank them with
> `lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")`.
>
> Write a short block of base-R code that independently verifies the *direction* of the
> coefficient without using `results()` again. It should take the gene with the largest absolute
> shrunken log2 fold change, pull DESeq2 median-of-ratios normalized counts with
> `counts(dds, normalized = TRUE)`, compute a naive log2 ratio of the mean normalized count in
> treated over the mean in control, print the control mean, the treated mean, the naive ratio,
> the unshrunken log2FC and the shrunken log2FC side by side, and `stopifnot()` that the naive
> ratio and the DESeq2 log2FC have the same sign. Add a small pseudocount so the ratio cannot
> divide by zero. Do not use any package other than DESeq2.

## 3. What AI generated

The block now present in `week5_deseq2_analysis.R` under
`# 7b. Independent verification of coefficient direction`: it computes `norm_counts`, selects
`top_gene` by `which.max(abs(res_df$log2FoldChange))`, indexes samples with
`colData(dds)$condition`, forms `log2((mean(treated) + 0.5) / (mean(control) + 0.5))`, prints the
comparison, and asserts sign agreement with `stopifnot()`.

## 4. Generated code was run locally

Yes. Run with `Rscript week5_deseq2_analysis.R` on R 4.3.3 / DESeq2 1.42.0 / apeglm 1.27.0.
Console output:

```
-- Direction check on Gene0008 --
  mean normalized count, control: 79.9
  mean normalized count, treated: 317.8
  naive log2 ratio (treated/control): 1.984
  DESeq2 unshrunken log2FC        : 2.018
  apeglm shrunken log2FC          : 1.933
```

## 5. Package functions and arguments I checked myself

| Item | How I checked | Result |
|---|---|---|
| `counts(dds, normalized = TRUE)` returns median-of-ratios normalized counts, not raw | `?counts` (DESeq2 method) | Confirmed; requires size factors, which `DESeq()` has already estimated |
| `lfcShrink(..., type = "apeglm")` needs `coef=`, not `contrast=` | `?lfcShrink` | Confirmed — apeglm only supports `coef`; this is why the script uses `contrast=` for `results()` and `coef=` for shrinkage |
| The coefficient name is not assumed | `resultsNames(dds)` printed in the script | Returned `"Intercept" "batch_B_vs_A" "batch_C_vs_A" "condition_treated_vs_control"`; the script `grep`s for the exact name and `stop()`s if absent |
| `results()` reports the direction I intended | `mcols(res)$description[2]` | `"log2 fold change (MLE): condition treated vs control"` |
| `vst(..., blind = FALSE)` is appropriate for design-aware QC | `?vst` | Confirmed |

## 6. Sample identity — independently verified

Not delegated to AI. In the script, before any modelling:

- `setequal()` plus `identical(colnames(counts_mat), rownames(coldata))` after explicit reordering,
  so a metadata file in a different order cannot silently mislabel samples;
- `duplicated()` checks on both column names and metadata row names;
- `all(counts_mat >= 0)` and `all(counts_mat == round(counts_mat))` to confirm the matrix is raw
  non-negative integers (no TPM/CPM/z-score);
- `identical(levels(coldata$condition), c("control", "treated"))` and
  `identical(levels(coldata$batch), c("A", "B", "C"))` to confirm both factors have the expected
  levels, plus `levels(coldata$condition)[1] == "control"` for the reference level;
- `qr(model.matrix(~ batch + condition))$rank` equals the number of columns (4 of 4), confirming a
  full-rank, non-confounded design;
- `table(batch, condition)` printed: 2 samples in each of the six cells, i.e. balanced.

I also checked the sample naming convention against the metadata by eye: `CA1/CA2` are control
batch A, `TA1/TA2` treated batch A, and so on for B and C — the names agree with the `condition`
and `batch` columns for all 12 samples.

## 7. Coefficient direction — independently verified

I did not take the AI block's `stopifnot()` on faith. I recomputed the top gene's ratio by hand
from the printed means:

```
log2((317.8 + 0.5) / (79.9 + 0.5)) = log2(318.3 / 80.4) = log2(3.959) = 1.985
```

which matches the script's `1.984` to rounding, and has the same positive sign as both the
unshrunken (`2.018`) and shrunken (`1.933`) DESeq2 estimates. So a positive value in
`week5_deseq2_results.csv` means **higher in treated**, and the `direction` column is labelled
accordingly.

As a second, independent check on direction I confirmed that the volcano plot's red points
("Up in treated") all lie at positive log2 fold change and the blue points at negative, which
would be visibly wrong had the contrast been inverted.

## 8. Errors and revisions — documented

**Error 1 (a starter-code defect that AI did not catch).** To be precise about attribution: the
failing line `vsd <- vst(dds, blind = FALSE)` is verbatim from the provided
`Week5_Homework_Starter.R` — it is not AI-generated. It is logged here because I asked the
assistant to sanity-check my script before I ran it and it raised no objection to this line, so
the error survived to runtime. On this dataset it failed:

```
Error in vst(dds, blind = FALSE) : less than 'nsub' rows,
  it is recommended to use varianceStabilizingTransformation directly
```

Cause, which I traced myself with `?vst`: `vst()` defaults to `nsub = 1000`, but only **989**
genes survive the `>= 10 counts in >= 3 samples` filter, so there were fewer rows than the
function wanted to subsample. **Revision:** `vst(dds, blind = FALSE, nsub = min(1000, nrow(dds)))`,
with a comment recording why. The script then ran to completion. This is a good illustration of
why code has to be executed rather than trusted, whether it comes from a starter template or from
an assistant — the call is syntactically fine and looks correct.

**Error 2 (my prompt's fault, caught before use).** My first, shorter prompt asked only for
"a check that the fold change sign is right", and the AI proposed comparing raw (unnormalized)
column means. I rejected that: library sizes here range from 112,619 to 151,748, a ~35% spread,
so raw means confound expression with sequencing depth. I rewrote the prompt to require
`counts(dds, normalized = TRUE)`, which is the version preserved in Section 2 above.

## 9. What AI did *not* do

- It did not write any part of the script other than Section 7b.
- It did not choose the design formula, the filtering rule, or the thresholds.
- It did not write the interpretation in `week5_interpretation.md`.
- It did not produce any of the numbers reported in that interpretation; those come from the
  script's own console output and from `week5_deseq2_results.csv`.
- It was not asked to judge whether any gene is biologically important.
