suppressPackageStartupMessages({
  library(DESeq2)
})
counts_file <- "C:/Users/Lenovo/workspace/week5_hw2_out/rnaseq_counts.tsv"
group_file  <- "C:/Users/Lenovo/workspace/week5_hw2_out/rnaseq_groups.tsv"
suppressPackageStartupMessages(library(data.table))
ct <- data.table::fread(counts_file, sep="\t", header=TRUE, check.names=FALSE)
cat("ct dim:", paste(dim(ct), collapse="x"), "  class:", class(ct)[1], "\n")
gene_ids <- ct$gene
ct[["gene"]] <- NULL
counts <- as.matrix(ct)
storage.mode(counts) <- "integer"
rownames(counts) <- gene_ids
cat("counts dim:", paste(dim(counts), collapse="x"), "\n")
meta <- read.table(group_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
cat("meta rows:", nrow(meta), "\n")
counts <- counts[, meta$sample]
finite_ok <- rowSums(is.finite(counts) & counts > 0) == ncol(counts)
cat("rows kept:", sum(finite_ok), "/", nrow(counts), "\n")
counts <- counts[finite_ok, ]
meta$group <- factor(meta$group)
# use sample IDs as rownames of colData (DESeqDataSetFromMatrix requires alignment)
rownames(meta) <- meta$sample
meta <- meta[colnames(counts), , drop=FALSE]
storage.mode(counts) <- "integer"
dds <- DESeqDataSetFromMatrix(countData=counts, colData=meta, design=~group)
dds <- DESeq(dds, quiet=TRUE)
cat("resultsNames:\n"); print(resultsNames(dds))

tests <- intersect(c("T4400+LIPUS","T3976","DMSO+LIPUS","T4400","T3976+LIPUS"),
                   levels(meta$group))
cat("tests:", paste(tests, collapse=","), "\n")
out_dir <- "C:/Users/Lenovo/workspace/week5_hw2_out"
for (t in tests){
  res <- results(dds, contrast=c("group", t, "DMSO"))
  df  <- as.data.frame(res)
  df$feature <- rownames(df)
  df$contrast <- paste0(t,"_vs_DMSO")
  fn <- file.path(out_dir, sprintf("diff_%s_vs_DMSO.csv", t))
  write.csv(df, fn, row.names=FALSE)
  cat("wrote:", fn, "  sig(padj<0.05 & |LFC|>=1) =",
      sum(df$padj<0.05 & abs(df$log2FoldChange)>=1, na.rm=TRUE), "\n")
}
cat("DONE\n")