# CircTest (dieterich-lab) on DCC's CircRNACount and LinearCount for one A-vs-B comparison.
#
# Filtering follows Circ.filter: at least filter_sample samples with filter_count circRNA reads,
# and circ / (circ + linear) of at least `percentage` in one group. Circ.filter itself is called
# when both groups have the same size (it splits columns into equal blocks); with unequal groups
# the same rule is applied per group here. Circ.test (beta-binomial, aod::betabin; null ~ 1 vs
# ~ group) is then called on each kept circRNA alone, so one failed fit only loses that circRNA;
# p-values are BH-adjusted over the tested circRNAs, as Circ.test adjusts over its input.
#
# Usage: Rscript circrna_circtest.R CircRNACount.tsv.gz LinearCount.tsv.gz n_a n_b
#                filter_count filter_sample percentage circtest_lib out.tsv

args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(args[8], .libPaths()))
suppressPackageStartupMessages(library(CircTest))

circ <- read.delim(gzfile(args[1]), check.names = FALSE, as.is = TRUE)
linear <- read.delim(gzfile(args[2]), check.names = FALSE, as.is = TRUE)
n_a <- as.integer(args[3]); n_b <- as.integer(args[4])
filter_count <- as.numeric(args[5]); filter_sample <- as.numeric(args[6]); percentage <- as.numeric(args[7])

if (nrow(circ) == 0) {
    writeLines(paste(c(names(circ)[1:3], "circ_reads", "ratio_a", "ratio_b", "p_value", "q_value", "status"),
                     collapse = "\t"), args[9])
    quit(save = "no")
}
if (nrow(circ) != nrow(linear) || any(circ[, 1:3] != linear[, 1:3]))
    stop("CircRNACount and LinearCount rows differ")
if (ncol(circ) - 3 != n_a + n_b || ncol(linear) - 3 != n_a + n_b)
    stop("sample columns do not match the comparison's replicates")
# DCC's CircRNACount and LinearCount both start with Chr, Start, End (CircTest's default description)
names(linear) <- names(circ)
description <- 1:3
group <- c(rep(1, n_a), rep(2, n_b))
counts <- as.matrix(circ[, -description, drop = FALSE])
host <- as.matrix(linear[, -description, drop = FALSE])

passes <- function(i) {
    if (sum(counts[i, ] >= filter_count) < filter_sample) return(FALSE)
    perc <- sapply(1:2, function(g) {
        c_sum <- sum(counts[i, group == g]); l_sum <- sum(host[i, group == g])
        c_sum / (c_sum + l_sum)
    })
    isTRUE(max(perc, na.rm = TRUE) >= percentage)
}
if (n_a == n_b) {
    kept_rows <- Circ.filter(circ = circ, linear = linear, Nreplicates = n_a, filter.sample = filter_sample,
                             filter.count = filter_count, percentage = percentage,
                             circle_description = description)
    kept <- rownames(circ) %in% rownames(kept_rows)
} else {
    kept <- vapply(seq_len(nrow(circ)), passes, logical(1))
}

ratio_mean <- function(i, g) {
    r <- counts[i, group == g] / (counts[i, group == g] + host[i, group == g])
    mean(r[is.finite(r)])
}
p_value <- rep(NA_real_, nrow(circ))
status <- ifelse(kept, "tested", "filtered")
for (i in which(kept)) {
    p <- tryCatch(Circ.test(Circ = circ[i, , drop = FALSE], Linear = linear[i, , drop = FALSE],
                            CircCoordinates = circ[i, description, drop = FALSE], group = group,
                            plotsig = FALSE, circle_description = description)$p.val,
                  error = function(e) { message("fit failed for row ", i, ": ", conditionMessage(e)); NA })
    if (length(p) != 1 || is.na(p)) status[i] <- "fit_failed" else p_value[i] <- p
}
q_value <- rep(NA_real_, nrow(circ))
tested <- which(status == "tested")
q_value[tested] <- p.adjust(p_value[tested], method = "BH")

out <- data.frame(circ[, description], circ_reads = rowSums(counts),
                  ratio_a = sapply(seq_len(nrow(circ)), ratio_mean, g = 1),
                  ratio_b = sapply(seq_len(nrow(circ)), ratio_mean, g = 2),
                  p_value = p_value, q_value = q_value, status = status, check.names = FALSE)
write.table(out, args[9], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message(sum(status == "tested"), " tested, ", sum(status == "filtered"), " filtered, ",
        sum(status == "fit_failed"), " fit failed")
