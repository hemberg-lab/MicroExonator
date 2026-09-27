#!/usr/bin/env Rscript
# Differential expression for one comparison, design ~ group (group_a vs group_b).
#
# Two routes, one script:
#   --route featurecounts  joined featureCounts matrix (gene_id, length, one
#                          column per biological replicate, technical runs summed)
#   --route tximport       joined Salmon gene counts (per replicate), TPM and
#                          effective length (per run; averaged per replicate here),
#                          through DESeqDataSetFromTximport
#
# When the preflight says inference is not supported, only descriptive output
# is written (size-factor normalised means per group); the status file gives
# the reasons and no test is run.

suppressPackageStartupMessages({
  library(DESeq2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list()
for (i in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[i])]] <- args[i + 1]
for (required in c("route", "preflight", "out")) {
  if (is.null(opt[[required]])) stop("missing --", required)
}

preflight <- fromJSON(opt$preflight)
replicates <- c(preflight$replicates$a, preflight$replicates$b)
group <- factor(c(rep("a", length(preflight$replicates$a)), rep("b", length(preflight$replicates$b))),
                levels = c("b", "a"))
col_data <- data.frame(group = group, row.names = replicates)

read_matrix <- function(path, key_columns) {
  table <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  values <- as.matrix(table[, -seq_len(key_columns), drop = FALSE])
  rownames(values) <- table[[1]]
  values
}

per_replicate_mean <- function(values) {
  # columns are runs; average the runs of each biological replicate
  collapse <- unlist(preflight$collapse)
  sapply(replicates, function(r) rowMeans(values[, names(collapse)[collapse == r], drop = FALSE]))
}

if (opt$route == "featurecounts") {
  counts <- read_matrix(opt$counts, 2)[, replicates, drop = FALSE]
  dds <- DESeqDataSetFromMatrix(round(counts), col_data, ~ group)
} else if (opt$route == "tximport") {
  counts <- read_matrix(opt$counts, 1)[, replicates, drop = FALSE]
  genes <- rownames(counts)
  abundance <- per_replicate_mean(read_matrix(opt$tpm, 1))[genes, , drop = FALSE]
  lengths <- per_replicate_mean(read_matrix(opt$length, 1))[genes, , drop = FALSE]
  txi <- list(abundance = abundance, counts = counts, length = lengths,
              countsFromAbundance = "no")
  dds <- DESeqDataSetFromTximport(txi, col_data, ~ group)
} else {
  stop("unknown route: ", opt$route)
}

dds <- estimateSizeFactors(dds)
normalised <- counts(dds, normalized = TRUE)
descriptive <- data.frame(gene_id = rownames(normalised),
                          mean_a = rowMeans(normalised[, group == "a", drop = FALSE]),
                          mean_b = rowMeans(normalised[, group == "b", drop = FALSE]))
write.table(descriptive, gzfile(paste0(opt$out, ".descriptive.tsv.gz")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Fixed outputs for the workflow: results (a comment line only when no test is
# run) and a one-line status.
results_path <- paste0(opt$out, ".results.tsv.gz")
if (isTRUE(preflight$inference_supported)) {
  dds <- DESeq(dds, quiet = TRUE)
  result <- results(dds, contrast = c("group", "a", "b"))
  result <- data.frame(gene_id = rownames(result), as.data.frame(result))
  write.table(result[order(result$padj), ], gzfile(results_path),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  writeLines("supported", paste0(opt$out, ".status.txt"))
} else {
  writeLines("# inference not supported", gzfile(results_path))
  writeLines(c("unsupported", unlist(preflight$reasons)), paste0(opt$out, ".status.txt"))
}
