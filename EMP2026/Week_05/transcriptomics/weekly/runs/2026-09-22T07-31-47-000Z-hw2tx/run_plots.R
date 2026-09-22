suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})
out_dir <- "C:/Users/Lenovo/workspace/week5_hw2_out"
counts_file <- file.path(out_dir, "rnaseq_counts.tsv")
group_file  <- file.path(out_dir, "rnaseq_groups.tsv")

suppressPackageStartupMessages(library(data.table))
ct <- data.table::fread(counts_file, sep="\t", header=TRUE, check.names=FALSE)
gene_ids <- ct$gene; ct[["gene"]] <- NULL
counts <- as.matrix(ct); storage.mode(counts) <- "integer"; rownames(counts) <- gene_ids
meta <- read.table(group_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
counts <- counts[, meta$sample]
finite_ok <- rowSums(is.finite(counts) & counts > 0) == ncol(counts)
counts <- counts[finite_ok, ]
rownames(meta) <- meta$sample
meta <- meta[colnames(counts), , drop=FALSE]
meta$group <- factor(meta$group)

# VST for visualization (sized speed advantage)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta, design=~group)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind=TRUE)

# PCA
pca <- plotPCA(vsd, intgroup="group", returnData=TRUE)
pct <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color=group, label=name)) +
  geom_point(size=3) + geom_text(vjust=-1, size=3) +
  scale_color_brewer(palette="Set2") +
  labs(title="PCA of RNA-seq samples (VST)",
       x=paste0("PC1 (", pct[1], "%)"), y=paste0("PC2 (", pct[2], "%)")) +
  theme_minimal(base_size=12) + theme(legend.position="bottom")
ggsave(file.path(out_dir, "pca.png"), p, width=8, height=5, units="in", bg="white")

# Heatmap of top 50 most variable
top_n_genes <- 50
rv <- apply(assay(vsd), 1, var)
W <- assay(vsd)[order(rv, decreasing=TRUE)[1:top_n_genes], ]
ann <- data.frame(group=meta$group); rownames(ann) <- meta$sample
png(file.path(out_dir, "heatmap_topvar50.png"), width=900, height=1100)
pheatmap(W, scale="row", annotation_col=ann, show_rownames=TRUE, fontsize_row=6,
         clustering_method="complete", color=colorRampPalette(c("navy","white","firebrick"))(100))
dev.off()

# Volcano for T4400 vs DMSO (largest effect)
res <- read.csv(file.path(out_dir, "diff_T4400_vs_DMSO.csv"))
res$sig <- ifelse(is.na(res$padj),"NS",
                  ifelse(res$padj<0.05 & abs(res$log2FoldChange)>=1, "sig", "ns"))
res$sign <- ifelse(res$sig=="sig", ifelse(res$log2FoldChange>0, "Up","Down"), "ns")
res <- res[!is.na(res$padj), ]
pv <- ggplot(res, aes(log2FoldChange, -log10(padj), color=sign)) +
  geom_point(alpha=0.7, size=1.5) +
  scale_color_manual(values=c("Up"="firebrick","Down"="steelblue","ns"="grey70")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  labs(title="Volcano: T4400 vs DMSO (DESeq2)",
       x="log2 fold change", y="-log10(adjusted p-value)") +
  theme_minimal(base_size=12) + theme(legend.position="bottom")
ggsave(file.path(out_dir, "volcano_T4400.png"), pv, width=7, height=5, units="in", bg="white")

# Top-20 boxplot for T4400 vs DMSO
top20 <- res[order(res$padj), ][1:20, c("feature","log2FoldChange","padj")]
vst_mat <- assay(vsd)[top20$feature, ]
df <- as.data.frame(t(vst_mat))
df$sample <- rownames(df)
df$group  <- meta[rownames(df), "group"]
df_long <- reshape2::melt(df, id.vars=c("sample","group"))
names(df_long) <- c("sample","group","gene","vst")
bp <- ggplot(df_long, aes(group, vst, fill=group)) +
  geom_boxplot(outlier.size=0.5) + facet_wrap(~gene, scales="free_y", ncol=5) +
  scale_fill_brewer(palette="Set2") +
  labs(title="Top-20 DEGs: T4400 vs DMSO (VST counts)",
       x=NULL, y="VST") +
  theme_minimal(base_size=10) + theme(legend.position="none",
                                       strip.text=element_text(size=6))
ggsave(file.path(out_dir, "top20_boxplot_T4400.png"), bp, width=12, height=8, units="in", bg="white")

# Save top DEG list
top_all <- do.call(rbind, lapply(c("T4400_vs_DMSO","T4400+LIPUS_vs_DMSO"), function(cn){
  df <- read.csv(file.path(out_dir, sprintf("diff_%s.csv", cn)))
  df <- df[order(df$padj), ]
  df <- df[!is.na(df$padj), ][1:30, ]
  df$contrast <- cn
  df[, c("feature","contrast","baseMean","log2FoldChange","pvalue","padj")]
}))
write.csv(top_all, file.path(out_dir, "top30_DEGs.csv"), row.names=FALSE)

cat("DONE\n")