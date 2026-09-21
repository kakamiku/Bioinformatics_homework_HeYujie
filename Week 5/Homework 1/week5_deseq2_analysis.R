# ================================================================
# Week 5 Homework 1
# From Verified Counts to an Interpretable DESeq2 Result
# Bulk RNA-seq differential expression with DESeq2
#
# Course : Bioinformatics: From Multi-Omics Data to Discovery
# Design : ~ batch + condition   (3 balanced batches A/B/C; control vs treated)
# ================================================================

# Required packages (install once):
# BiocManager::install(c("DESeq2", "apeglm"))
# install.packages(c("tidyverse", "ggrepel"))

# Set working directory to the folder that holds this script.
# (Interactive RStudio convenience; skipped when run with Rscript.)
if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
})

set.seed(20260921)   # reproducibility for anything stochastic

# ----------------------------
# 1. Project paths
# ----------------------------
count_file    <- "Week5_Homework_Count_Matrix.csv"
metadata_file <- "Week5_Homework_Sample_Metadata.csv"

# All eight submission files are written flat, into the same folder as this
# script, under exactly the filenames the assignment asks for.

# Analysis parameters, declared once and reported in every output.
ALPHA      <- 0.05   # FDR (adjusted p value) threshold
LFC_CUT    <- 1      # |log2 fold change| effect-size threshold
MIN_COUNT  <- 10     # pre-filter: minimum count ...
MIN_SAMPLE <- 3      # ... in at least this many samples

# ----------------------------
# 2. Import (inputs are read only; nothing is written back to them)
# ----------------------------
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)

coldata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)

cat("Imported count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

# ----------------------------
# 3. Mandatory validation
# ----------------------------
counts_mat <- as.matrix(counts)

# 3a. Values are non-negative integers (no TPM/CPM/z-score transformation).
stopifnot(is.numeric(counts_mat))
stopifnot(!any(is.na(counts_mat)))
stopifnot(all(counts_mat >= 0))
stopifnot(all(counts_mat == round(counts_mat)))
storage.mode(counts_mat) <- "integer"

# 3b. Sample identity: same number, same names, same order, no duplicates.
stopifnot(ncol(counts_mat) == nrow(coldata))
stopifnot(!any(duplicated(colnames(counts_mat))))
stopifnot(!any(duplicated(rownames(coldata))))
stopifnot(setequal(colnames(counts_mat), rownames(coldata)))
# Reorder metadata to the matrix column order, then require exact identity.
coldata <- coldata[colnames(counts_mat), , drop = FALSE]
stopifnot(identical(colnames(counts_mat), rownames(coldata)))

# 3c. Factors with the intended reference level.
coldata$condition <- relevel(factor(coldata$condition), ref = "control")
coldata$batch     <- factor(coldata$batch)
stopifnot(is.factor(coldata$condition), is.factor(coldata$batch))
stopifnot(levels(coldata$condition)[1] == "control")   # control is the reference
stopifnot(identical(levels(coldata$condition), c("control", "treated")))
stopifnot(identical(levels(coldata$batch), c("A", "B", "C")))      # expected batch levels

cat("\n-- Design balance (batch x condition) --\n")
print(table(coldata$batch, coldata$condition))
cat("\n-- Library sizes (total counts per sample) --\n")
print(summary(colSums(counts_mat)))

# 3d. The model matrix must be full rank for ~ batch + condition.
mm <- model.matrix(~ batch + condition, data = coldata)
stopifnot(qr(mm)$rank == ncol(mm))
cat("\nModel matrix is full rank:", qr(mm)$rank, "of", ncol(mm), "columns\n")

# ----------------------------
# 4. Construct DESeq2 object
# ----------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = coldata,
  design    = ~ batch + condition
)

# Why batch is in the design:
#   Samples were processed in three batches (A, B, C), each containing one
#   control/treated pair per group, so batch is balanced with respect to
#   condition. Technical batch effects shift a gene's expression in every
#   sample of that batch. Leaving batch out would push that systematic shift
#   into the residual (dispersion) term, inflating within-group variance and
#   costing power; including it as a blocking factor lets DESeq2 estimate and
#   remove the batch offsets, so the condition coefficient reports the
#   treatment effect *after* adjusting for batch. Because the design is
#   balanced, batch and condition are not confounded and the model is full
#   rank (checked above).

# ----------------------------
# 5. Pre-filter
# ----------------------------
# Filtering rule (reported): keep a gene if it has at least MIN_COUNT counts
# in at least MIN_SAMPLE samples. MIN_SAMPLE = 3 is smaller than the smallest
# group of interest (6 per condition), so a gene expressed in only one
# condition is retained.
keep <- rowSums(counts(dds) >= MIN_COUNT) >= MIN_SAMPLE
cat("\nFilter rule: >=", MIN_COUNT, "counts in >=", MIN_SAMPLE, "samples\n")
cat("Genes before filtering:", nrow(dds), "\n")
dds <- dds[keep, ]
cat("Genes after filtering :", nrow(dds), "\n")

# ----------------------------
# 6. Fit model
# ----------------------------
dds <- DESeq(dds)

coef_names <- resultsNames(dds)
cat("\n-- resultsNames(dds) --\n")
print(coef_names)

# Do not assume the coefficient name: confirm it against resultsNames(dds).
target_coef <- grep("^condition_treated_vs_control$", coef_names, value = TRUE)
if (length(target_coef) != 1) {
  stop("Expected coefficient was not found. Inspect resultsNames(dds) and update target_coef.")
}
cat("Confirmed shrinkage coefficient:", target_coef, "\n")

# ----------------------------
# 7. Extract and shrink results
# ----------------------------
res <- results(
  dds,
  contrast = c("condition", "treated", "control"),
  alpha    = ALPHA
)
cat("\n-- results() metadata (confirms contrast direction) --\n")
print(mcols(res)$description[2])
cat("\n-- summary(res), unshrunken --\n")
summary(res)

res_shrunk <- lfcShrink(
  dds,
  coef = target_coef,
  type = "apeglm"
)

# Sanity check: shrinkage changes the magnitude of the LFC, never the contrast.
stopifnot(identical(rownames(res), rownames(res_shrunk)))

res_df <- as.data.frame(res_shrunk) |>
  rownames_to_column("gene_id") |>
  mutate(
    log2FoldChange_unshrunk = res$log2FoldChange,   # kept for transparency
    stat_wald               = res$stat,
    significant = !is.na(padj) & padj < ALPHA & abs(log2FoldChange) >= LFC_CUT,
    direction = case_when(
      significant & log2FoldChange > 0 ~ "Up in treated",
      significant & log2FoldChange < 0 ~ "Down in treated",
      TRUE                             ~ "Not significant"
    )
  ) |>
  arrange(padj, desc(abs(log2FoldChange)))

# The COMPLETE table is exported: every gene that survived filtering,
# significant or not.
#
# Column dictionary (note the mix of shrunken and unshrunken quantities, which
# is standard DESeq2 practice):
#   baseMean                - mean of normalized counts across all samples
#   log2FoldChange, lfcSE   - apeglm SHRUNKEN posterior estimate and its SE
#   pvalue, padj, stat_wald - from the UNSHRUNKEN Wald test (apeglm does not
#                             recompute p values; shrinkage affects effect
#                             size only)
#   log2FoldChange_unshrunk - the raw MLE estimate, kept for comparison
#   significant, direction  - padj < ALPHA AND |SHRUNKEN log2FC| >= LFC_CUT
write.csv(res_df, "week5_deseq2_results.csv", row.names = FALSE)
cat("\nExported", nrow(res_df), "genes (all filtered genes, not only significant ones)\n")

cat("\nThresholds: padj <", ALPHA, "AND |log2FoldChange| >=", LFC_CUT,
    "(applied to the apeglm-SHRUNKEN log2FC)\n")
cat("  for reference - padj <", ALPHA, "alone            :",
    sum(!is.na(res_df$padj) & res_df$padj < ALPHA), "genes\n")
cat("  for reference - padj + |UNSHRUNKEN log2FC| >= 1 :",
    sum(!is.na(res_df$padj) & res_df$padj < ALPHA &
          abs(res_df$log2FoldChange_unshrunk) >= LFC_CUT), "genes\n")
cat("Significant genes:", sum(res_df$significant), "\n")
print(table(res_df$direction))

# ----------------------------
# 7b. Independent verification of coefficient direction
#     (see week5_AI_verification_log.md - AI drafted this block,
#      the numbers below were re-derived by hand from normalized counts)
# ----------------------------
norm_counts <- counts(dds, normalized = TRUE)
top_gene <- res_df$gene_id[which.max(abs(res_df$log2FoldChange))[1]]
ctrl_idx <- which(colData(dds)$condition == "control")
trt_idx  <- which(colData(dds)$condition == "treated")
manual_lfc <- log2((mean(norm_counts[top_gene, trt_idx]) + 0.5) /
                   (mean(norm_counts[top_gene, ctrl_idx]) + 0.5))
cat("\n-- Direction check on", top_gene, "--\n")
cat("  mean normalized count, control:", round(mean(norm_counts[top_gene, ctrl_idx]), 1), "\n")
cat("  mean normalized count, treated:", round(mean(norm_counts[top_gene, trt_idx]), 1), "\n")
cat("  naive log2 ratio (treated/control):",
    round(manual_lfc, 3), "\n")
cat("  DESeq2 unshrunken log2FC        :",
    round(res_df$log2FoldChange_unshrunk[res_df$gene_id == top_gene], 3), "\n")
cat("  apeglm shrunken log2FC          :",
    round(res_df$log2FoldChange[res_df$gene_id == top_gene], 3), "\n")
# A positive sign in all three means "higher in treated", i.e. the reported
# contrast really is treated versus control and not the reverse.
stopifnot(sign(manual_lfc) ==
            sign(res_df$log2FoldChange_unshrunk[res_df$gene_id == top_gene]))

# ----------------------------
# 8. PCA (on variance-stabilized data, not raw counts)
# ----------------------------
# nsub must not exceed the number of genes left after filtering (989 here),
# otherwise vst() errors; the default is 1000.
vsd <- vst(dds, blind = FALSE, nsub = min(1000, nrow(dds)))

pca_df <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = condition, shape = batch, label = name)
) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = c("control" = "#2F6DB3", "treated" = "#C0392B")) +
  labs(
    title    = "Week 5 RNA-seq PCA (VST, blind = FALSE)",
    subtitle = paste0(nrow(dds), " genes after filtering; colour = condition, shape = batch"),
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_bw(base_size = 12)

ggsave("week5_pca.png", p_pca, width = 7, height = 5, dpi = 300)

# ----------------------------
# 9. Volcano plot
# ----------------------------
n_up   <- sum(res_df$direction == "Up in treated")
n_down <- sum(res_df$direction == "Down in treated")

label_df <- res_df |>
  filter(significant) |>
  slice_head(n = 10)

p_volcano <- ggplot(
  res_df |> mutate(neg_log10_padj = -log10(pmax(padj, 1e-300))),
  aes(x = log2FoldChange, y = neg_log10_padj, color = direction)
) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed") +
  geom_hline(yintercept = -log10(ALPHA), linetype = "dashed") +
  geom_text_repel(
    data = label_df |> mutate(neg_log10_padj = -log10(pmax(padj, 1e-300))),
    aes(label = gene_id), size = 2.8, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Up in treated"   = "#C0392B",
    "Down in treated" = "#2F6DB3",
    "Not significant" = "grey70"
  )) +
  labs(
    title    = "Differential expression: treated versus control",
    subtitle = paste0(
      "Dashed lines: padj = ", ALPHA, " and |log2FC| = ", LFC_CUT,
      "  |  ", n_up, " up, ", n_down, " down of ", nrow(res_df), " tested"
    ),
    x     = "Shrunken log2 fold change (apeglm), treated / control",
    y     = "-log10 adjusted p value (BH)",
    color = NULL
  ) +
  theme_bw(base_size = 12)

ggsave("week5_de_plot.png", p_volcano, width = 7, height = 5, dpi = 300)

# ----------------------------
# 10. Save reproducibility files
# ----------------------------
saveRDS(dds, "week5_deseq2_object.rds")
capture.output(sessionInfo(), file = "session_info.txt")

cat("\nDone. All submission files written to the current directory.\n")

# ----------------------------
# 11. Interpretation
# ----------------------------
# See week5_interpretation.md (100-150 words) and week5_AI_verification_log.md.
