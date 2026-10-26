# Replicated BSJ abundance, using edgeR's quasi-likelihood workflow.
# This is not native CIRI_DE_replicate: normalization is TMMwsp on BSJ counts.
# Every outcome is a status, never a crash: insufficient data (too few replicates,
# events or reads after filterByExpr) and fit failures are reported with their reason,
# and R warnings are kept in status.tsv.
suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) %in% c(1, 2))
folder <- args[1]
min_events <- if (length(args) == 2) as.integer(args[2]) else 10L
counts <- read.delim(gzfile(file.path(folder, 'counts.tsv.gz')), check.names=FALSE)
samples <- read.delim(file.path(folder, 'samples.tsv'), stringsAsFactors=FALSE)
stopifnot(identical(names(counts)[-1], samples$sample), !anyDuplicated(samples$sample))
group <- factor(samples$group, levels=c('b', 'a'))
stopifnot(!anyNA(group))
status <- 'insufficient_data'
reason <- ''
warned <- character()
answer <- data.frame(circ_id=character(), logFC=numeric(), logCPM=numeric(), F=numeric(), PValue=numeric(), FDR=numeric())
x <- as.matrix(counts[,-1,drop=FALSE])
rownames(x) <- counts[[1]]
stopifnot(!anyDuplicated(rownames(x)), all(x[!is.na(x)] >= 0), all(x[!is.na(x)] == floor(x[!is.na(x)])))
x <- x[complete.cases(x),,drop=FALSE]
tested <- 0L
if (!all(table(group) >= 2)) {
    reason <- 'fewer than two biological replicates on a side'
} else if (nrow(x) < min_events) {
    reason <- sprintf('%d measured events, fewer than %d', nrow(x), min_events)
} else if (!all(colSums(x) > 0)) {
    reason <- 'a replicate has no BSJ read'
} else {
    design <- model.matrix(~group)
    y <- DGEList(x)
    keep <- filterByExpr(y, design=design)
    y <- y[keep,,keep.lib.sizes=FALSE]
    tested <- nrow(y)
    if (nrow(y) < min_events) {
        reason <- sprintf('%d events pass filterByExpr, fewer than %d', nrow(y), min_events)
    } else if (!all(colSums(y$counts) > 0)) {
        reason <- 'a replicate has no BSJ read among the filtered events'
    } else {
        fitted <- withCallingHandlers(tryCatch({
            y <- calcNormFactors(y, method='TMMwsp')
            y <- estimateDisp(y, design)
            fit <- glmQLFit(y, design, robust=TRUE)
            list(y=y, table=topTags(glmQLFTest(fit, coef=2), n=Inf, sort.by='none')$table)
        }, error=function(e) conditionMessage(e)),
        warning=function(w) {warned <<- c(warned, conditionMessage(w)); invokeRestart('muffleWarning')})
        if (is.character(fitted)) {
            status <- 'fit_failed'
            reason <- gsub('[\t\n]', ' ', fitted)
        } else {
            answer <- data.frame(circ_id=rownames(fitted$table), fitted$table, row.names=NULL)
            status <- 'tested'
            write.table(fitted$y$samples, file.path(folder, 'normalization.tsv'), sep='\t', quote=FALSE, col.names=NA)
        }
    }
}
con <- gzfile(file.path(folder, 'abundance.tsv.gz'), 'wt')
write.table(answer, con, sep='\t', quote=FALSE, row.names=FALSE)
close(con)
lines <- c(paste('status', status, sep='\t'), paste('reason', reason, sep='\t'),
           paste('events_measured', nrow(x), sep='\t'), paste('events_tested', tested, sep='\t'),
           paste('min_events', min_events, sep='\t'),
           paste('edgeR', packageVersion('edgeR'), sep='\t'),
           'effect\tlog2 abundance A/B', 'normalization\tTMMwsp on retained BSJ counts',
           'multiplicity\tBH within tool/comparison tested universe')
if (length(warned)) lines <- c(lines, paste('warning', gsub('[\t\n]', ' ', unique(warned)), sep='\t'))
writeLines(lines, file.path(folder, 'status.tsv'))
