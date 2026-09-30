#!/usr/bin/env Rscript
# QAPA site usage with DEXSeq: each poly(A) site is a "feature" of its gene;
# the interaction condition:feature tests a change in the site's share of the
# gene's reads (usage), not a change in its abundance. Counts are the
# per-replicate sums from run_qapa.py pau (technical runs already summed).
# Only genes with at least two sites are tested. Direction: A over B.
suppressPackageStartupMessages({
  library(DEXSeq)
  library(BiocParallel)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: qapa_dexseq.R site_counts.tsv.gz samples.tsv out.tsv.gz threads")
counts <- read.delim(gzfile(args[1]), check.names = FALSE)
samples <- read.delim(args[2], stringsAsFactors = FALSE)
threads <- as.integer(args[4])
multi <- names(which(table(counts$gene_id) >= 2))
if (min(table(samples$condition)) < 2) multi <- character(0)   # no test below 2 per side
counts <- counts[counts$gene_id %in% multi, ]
out <- gzfile(args[3], "w")
if (nrow(counts) == 0) {
  writeLines("site_id\tgene_id\tlog2fold_a_b\tpvalue\tpadj", out)
  close(out)
  quit(save = "no")
}
matrix <- round(as.matrix(counts[, samples$sample, drop = FALSE]))
storage.mode(matrix) <- "integer"
design <- data.frame(row.names = samples$sample,
                     condition = factor(samples$condition, levels = c("b", "a")))
dxd <- DEXSeqDataSet(matrix, sampleData = design, design = ~ sample + exon + condition:exon,
                     featureID = counts$site_id, groupID = counts$gene_id)
parallel <- if (threads > 1) MulticoreParam(threads) else SerialParam()
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = parallel)
dxd <- testForDEU(dxd, BPPARAM = parallel)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM = parallel)
result <- DEXSeqResults(dxd)
fold <- grep("^log2fold_", colnames(result), value = TRUE)[1]
table <- data.frame(site_id = result$featureID, gene_id = result$groupID,
                    log2fold_a_b = if (is.na(fold)) NA else result[[fold]],
                    pvalue = result$pvalue, padj = result$padj)
write.table(table, out, sep = "\t", quote = FALSE, row.names = FALSE)
close(out)
