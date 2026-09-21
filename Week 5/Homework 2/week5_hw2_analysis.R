# =============================================================
# Week 5 Homework 2 - PCA + Boxplot + DEG table
# Built on the EasyMultiProfiler-Web session 8qqwdkILKa26siMSZ4lnW3np
# =============================================================

set.seed(20260921)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(MultiAssayExperiment)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(matrixStats)
})

repo_root  <- normalizePath("C:/Users/Lenovo/EasyMultiProfiler-Web", winslash = "/", mustWork = FALSE)
local_libs <- file.path(repo_root, ".local_run", "R_libs")
if (dir.exists(local_libs)) .libPaths(c(local_libs, .libPaths()))

library(EasyMultiProfiler)

session_dir <- "C:/Users/Lenovo/EasyMultiProfiler-Web/.local_run/data/sessions/8qqwdkILKa26siMSZ4lnW3np"

# Use the 24-sample raw empt, not the subsetted empt from the last diff cache
raw_empt <- readRDS(file.path(session_dir, "raw_empt_rna_seq_bmmsc.rds"))
cat("Loaded raw_empt:", class(raw_empt)[1], "with", nrow(raw_empt), "genes x",
    ncol(raw_empt), "samples\n")
cd <- as.data.frame(SummarizedExperiment::colData(raw_empt))
cat("Groups:\n")
print(table(cd$Group))

assay_data <- SummarizedExperiment::assays(raw_empt)[[1]]
cat("Assay range:", range(assay_data), "\n")

# ----------------- Run DESeq2 directly on full 24-sample dataset -----------------
dds <- DESeqDataSetFromMatrix(
  countData = assay_data,
  colData   = cd,
  design    = ~ Group
)
keep <- rowSums(counts(dds) >= 10) >= 3
cat("Filter: >=10 counts in >=3 samples; kept", sum(keep), "of", nrow(dds), "genes\n")
dds <- dds[keep, ]
dds <- DESeq(dds)
cat("Coefficients:\n")
print(resultsNames(dds))

# T4400 vs DMSO - this is the contrast that produced the biggest signal
res <- results(dds, contrast = c("Group", "T4400", "DMSO"), alpha = 0.05)
cat("T4400 vs DMSO - summary:\n")
summary(res)

# VST on the full 24-sample dds for the PCA
vsd <- vst(dds, blind = FALSE, nsub = min(1000, nrow(dds)))

# Save DEG table (full, ordered)
res_df <- as.data.frame(res)
res_df$feature <- rownames(res_df)
res_df <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange)), ]
write.csv(res_df,
          "C:/Users/Lenovo/workspace/week5_hw2_deseq2_results.csv",
          row.names = FALSE)
cat("Saved week5_hw2_deseq2_results.csv\n")

sig_mask <- !is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1
cat("Significant (padj<0.05 & |log2FC|>=1):", sum(sig_mask), "\n")

# ----------------- PCA -----------------
vsd_mat <- SummarizedExperiment::assays(vsd)[[1]]
rv <- matrixStats::rowVars(vsd_mat)
top_var <- order(rv, decreasing = TRUE)[seq_len(min(1000, length(rv)))]
mat_top <- vsd_mat[top_var, ]
mat_top <- mat_top - rowMeans(mat_top)
pca <- prcomp(t(mat_top), center = FALSE, scale. = FALSE)
pv <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

pca_df <- data.frame(
  sample = colnames(vsd_mat),
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2],
  Group  = cd$Group,
  stringsAsFactors = FALSE
)

# Quick one-way ANOVA on PC1 by group
fit_pc1 <- summary(aov(PC1 ~ Group, data = pca_df))
cat("PC1 ~ Group ANOVA p:", fit_pc1[[1]]$`Pr(>F)`[1], "\n")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = sample)) +
  geom_point(size = 3.6) +
  geom_text_repel(size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title    = "Week 5 Homework 2 - RNA-seq PCA",
    subtitle = paste0("Top ", length(top_var), " variable genes, VST-transformed counts"),
    x = paste0("PC1: ", round(pv[1], 1), "% variance"),
    y = paste0("PC2: ", round(pv[2], 1), "% variance")
  ) +
  theme_bw(base_size = 12)

ggsave("C:/Users/Lenovo/workspace/week5_hw2_pca.png",
       p_pca, width = 8, height = 6, dpi = 200)
cat("Saved week5_hw2_pca.png\n")

# ----------------- Boxplot of top 10 DEGs -----------------
sig_df <- res_df[sig_mask, ]
sig_df <- sig_df[order(-abs(sig_df$log2FoldChange)), ]
top10 <- head(sig_df$feature, 10)
cat("Top 10 DEGs by |log2FC|:\n")
print(top10)

# Expression values for boxplot
vsd_top <- as.data.frame(t(vsd_mat[top10, , drop = FALSE]))
vsd_top$sample <- rownames(vsd_top)
vsd_top$Group  <- cd$Group
plot_df <- reshape2::melt(vsd_top, id.vars = c("sample", "Group"),
                           variable.name = "gene", value.name = "expr")
plot_df$gene <- factor(plot_df$gene, levels = top10)

p_box <- ggplot(plot_df, aes(x = Group, y = expr, fill = Group)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.7) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Week 5 Homework 2 - Top 10 DEGs (T4400 vs DMSO)",
    subtitle = "VST-transformed expression across all 6 conditions (4 reps each)",
    x = NULL, y = "VST expression"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

ggsave("C:/Users/Lenovo/workspace/week5_hw2_topdeg_box.png",
       p_box, width = 13, height = 7, dpi = 200)
cat("Saved week5_hw2_topdeg_box.png\n")

cat("\nDone.\n")
