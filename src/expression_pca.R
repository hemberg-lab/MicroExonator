#!/usr/bin/env Rscript
# Cached comparison QC only. Never removes samples or changes DE results.
suppressPackageStartupMessages({
  library(DESeq2)
  library(jsonlite)
  library(ggplot2)
  library(ggrepel)
})

read_pca_matrix <- function(path, keys) {
  tab <- read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
  if (anyDuplicated(tab[[1]]) || anyDuplicated(names(tab))) stop("Duplicate gene/sample identifiers")
  values <- as.matrix(tab[, -seq_len(keys), drop = FALSE])
  storage.mode(values) <- "double"
  rownames(values) <- tab[[1]]
  if (any(!is.finite(values)) || any(values < 0)) stop("Nonfinite or negative matrix values")
  values
}

group_distance_qc <- function(expression, groups, min_group = 5L) {
  score <- rep(NA_real_, length(groups))
  threshold <- rep(NA_real_, length(groups))
  status <- rep("Not assessed", length(groups))
  reason <- rep("Too few biological replicates", length(groups))
  for (g in unique(groups)) {
    idx <- which(groups == g)
    if (length(idx) < min_group) next
    distances <- as.matrix(dist(t(expression[, idx, drop = FALSE]))) / sqrt(nrow(expression))
    diag(distances) <- NA_real_
    values <- apply(distances, 1, median, na.rm = TRUE)
    score[idx] <- values
    spread <- mad(values)
    if (!is.finite(spread) || spread <= sqrt(.Machine$double.eps) * max(1, median(values))) {
      reason[idx] <- "Zero or near-zero distance dispersion"
      next
    }
    cutoff <- median(values) + 3 * spread
    threshold[idx] <- cutoff
    status[idx] <- ifelse(values > cutoff, "Review", "Unflagged")
    reason[idx] <- "Exploratory within-group median distance > median + 3 scaled MAD"
  }
  data.frame(distance_score = score, distance_threshold = threshold,
             qc_status = status, qc_reason = reason)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) %% 2) stop("Arguments must be --key value pairs")
  opt <- list()
  for (i in seq(1, length(args), by = 2)) opt[[sub("^--", "", args[i])]] <- args[i+1]
  for (key in c("route", "preflight", "counts", "out")) if (is.null(opt[[key]])) stop("Missing --", key)
  min_group <- if (is.null(opt[["min-group"]])) 5L else as.integer(opt[["min-group"]])
  if (is.na(min_group) || min_group < 5) stop("--min-group must be at least 5")
  pf <- fromJSON(opt$preflight)
  samples <- c(pf$replicates$a, pf$replicates$b)
  if (length(samples) < 3 || anyDuplicated(samples)) stop("PCA needs at least three distinct biological replicates")
  groups <- c(rep(pf$group_a, length(pf$replicates$a)), rep(pf$group_b, length(pf$replicates$b)))
  if (length(groups) != length(samples) || pf$group_a == pf$group_b) stop("Invalid comparison group labels")
  counts <- read_pca_matrix(opt$counts, if (opt$route == "featurecounts") 2L else 1L)
  missing <- setdiff(samples, colnames(counts))
  if (length(missing)) stop("Missing biological replicate columns: ", paste(missing, collapse = ", "))
  counts <- counts[, samples, drop = FALSE]
  raw_counts <- counts
  keep <- rowSums(counts >= 10) >= 2
  counts <- counts[keep, , drop = FALSE]
  if (nrow(counts) < 2 || any(colSums(counts) == 0)) stop("Insufficient expressed genes or empty sample")
  metadata <- data.frame(group = factor(groups), row.names = samples)
  if (opt$route == "featurecounts") {
    dds <- DESeqDataSetFromMatrix(round(counts), metadata, ~ group)
  } else if (opt$route == "tximport") {
    if (is.null(opt$tpm) || is.null(opt$length)) stop("Salmon PCA requires TPM and effective length caches")
    collapse <- unlist(pf$collapse)
    replicate_mean <- function(path) {
      values <- read_pca_matrix(path, 1L)
      if (!all(rownames(counts) %in% rownames(values))) stop("Salmon gene identifiers do not match")
      output <- vapply(samples, function(sample) {
        runs <- names(collapse)[collapse == sample]
        if (!length(runs) || !all(runs %in% colnames(values))) stop("Missing Salmon run columns for ", sample)
        rowMeans(values[rownames(counts), runs, drop = FALSE])
      }, numeric(nrow(counts)))
      rownames(output) <- rownames(counts)
      output
    }
    txi <- list(counts = counts, abundance = replicate_mean(opt$tpm),
                length = replicate_mean(opt$length), countsFromAbundance = "no")
    if (any(txi$length <= 0)) stop("Nonpositive Salmon effective lengths")
    dds <- DESeqDataSetFromTximport(txi, metadata, ~ group)
  } else stop("Unknown expression route: ", opt$route)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  # Full VST works on small gene panels too (vst() requires many expressed genes).
  transformed <- assay(varianceStabilizingTransformation(dds, blind = TRUE, fitType = "mean"))
  variable <- apply(transformed, 1, var)
  selected <- head(order(variable, decreasing = TRUE), 500L)
  selected <- selected[variable[selected] > 0]
  if (length(selected) < 2) stop("Insufficient variable genes for PCA")
  pca <- prcomp(t(transformed[selected, , drop = FALSE]), center = TRUE, scale. = FALSE)
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  scores <- data.frame(sample_id = samples, group = groups,
                       PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                       PC1_variance_pct = explained[1], PC2_variance_pct = explained[2],
                       group_distance_qc(transformed, groups, min_group), row.names = NULL)
  qc_report <- identical(opt[["qc-report"]], "true")
  if (qc_report) {
    scores$qc_status <- "Not classified"
    scores$distance_threshold <- NA_real_
    scores$qc_reason <- "Descriptive QC only; no automatic sample exclusions"
  }
  scores$qc_status <- factor(scores$qc_status, levels = c("Unflagged", "Review", "Not assessed", "Not classified"))
  scores$group <- factor(scores$group, levels = c(pf$group_a, pf$group_b))
  pca_plot <- ggplot(scores, aes(PC1, PC2, colour = group, shape = qc_status)) +
    geom_point(size = 3.2) +
    geom_text_repel(aes(label = sample_id), seed = 42, max.overlaps = Inf,
                    max.iter = 10000, max.time = 2, box.padding = 0.5,
                    min.segment.length = 0, show.legend = FALSE, size = 3.2) +
    scale_colour_manual(values = setNames(c("#0072B2", "#D55E00"), c(pf$group_a, pf$group_b))) +
    scale_shape_manual(values = c("Unflagged" = 16, "Review" = 17, "Not assessed" = 1)) +
    labs(title = if (is.null(pf$comparison_id)) "Expression PCA" else pf$comparison_id,
         subtitle = paste(opt$route, "|", length(samples), "biological replicates"),
         x = sprintf("PC1 (%.1f%%)", explained[1]), y = sprintf("PC2 (%.1f%%)", explained[2]),
         colour = "Comparison group", shape = "QC status",
         caption = "VST; up to 500 most variable genes. Review flags are exploratory, not exclusions.") +
    theme_classic(base_size = 12) + theme(legend.position = "bottom")
  dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
  if (qc_report) {
    # The standalone report opts out of the old exploratory binary flags.
    pca_plot$mapping$shape <- NULL
    pca_plot <- pca_plot + guides(shape = "none") +
      labs(caption = "VST; up to 500 most variable genes. Descriptive QC only; no automatic exclusions.")
    script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
    source(file.path(dirname(script_path), "gene_expression_qc.R"))
    gene_expression_qc_report(opt, pf, dds, raw_counts, transformed, pca_plot)
  } else {
    ggsave(paste0(opt$out, ".pdf"), plot = pca_plot, width = 9, height = 6, units = "in")
  }
  write.table(scores, gzfile(paste0(opt$out, ".scores.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
}

if (sys.nframe() == 0) main()
