# Multi-page descriptive QC from retained comparison matrices only.
gene_expression_qc_report <- function(opt, pf, dds, raw, transformed, original_pca) {
  # Grid layout may consult a graphics device before ggsave opens the final PDF.
  # A null device prevents an unintended Rplots.pdf in the run directory.
  pdf(file=NULL)
  on.exit(dev.off(),add=TRUE)
  suppressPackageStartupMessages(library(gridExtra))
  samples <- colnames(raw)
  groups <- as.character(colData(dds)$group)
  palette <- setNames(c("#0072B2", "#D55E00"), c(pf$group_a, pf$group_b))
  ranking <- order(apply(transformed, 1, var), decreasing = TRUE)
  ranking <- ranking[apply(transformed[ranking, , drop=FALSE], 1, var) > 0]
  # Gene selection ignores group labels and DE results.
  views <- list()
  for (n in unique(c(min(500,length(ranking)), min(2000,length(ranking)), length(ranking)))) {
    fit <- prcomp(t(transformed[head(ranking,n), , drop=FALSE]), scale.=FALSE)
    pct <- 100*fit$sdev^2/sum(fit$sdev^2)
    axes <- list(c(1,2))
    if (n == min(500,length(ranking)) && ncol(fit$x) >= 3) axes <- c(axes,list(c(2,3)))
    for (a in axes) {
      name <- sprintf("%s genes | PC%d %.1f%% / PC%d %.1f%%",
                      if(n==length(ranking)) "All variable filtered" else n,a[1],pct[a[1]],a[2],pct[a[2]])
      views[[length(views)+1]] <- data.frame(sample=samples,group=groups,
                                            x=fit$x[,a[1]],y=fit$x[,a[2]],panel=name)
    }
  }
  pc <- do.call(rbind,views)
  pc$panel <- factor(pc$panel,levels=unique(pc$panel))
  pca_views <- ggplot(pc,aes(x,y,colour=group,label=sample))+geom_point(size=2.5)+
    geom_text_repel(seed=42,max.overlaps=Inf,size=2.7,box.padding=.4,max.time=2)+
    facet_wrap(~panel,scales="free",ncol=2)+scale_colour_manual(values=palette)+theme_classic()+
    labs(title="PCA robustness and additional axes",x="PC score",y="PC score",colour="Group",
         caption="Same biological replicates; no batch removal. Separation is not proof of technical failure.")
  distance <- as.matrix(dist(t(transformed)))/sqrt(nrow(transformed))
  correlation <- cor(transformed[head(ranking,500), , drop=FALSE],method="pearson")
  ordered <- samples[hclust(as.dist(distance),method="average")$order]
  heatmap_plot <- function(values,title,legend,colours,limits=NULL) {
    frame <- expand.grid(x=seq_along(ordered),y=seq_along(ordered))
    frame$value <- as.vector(values[ordered,ordered])
    frame$y <- length(ordered)+1-frame$y
    annotation <- data.frame(x=seq_along(ordered),y=length(ordered)+1,
                             group=groups[match(ordered,samples)])
    ggplot(frame,aes(x,y,fill=value))+geom_tile()+
      geom_point(data=annotation,aes(x,y,colour=group),inherit.aes=FALSE,shape=15,size=4)+
      scale_fill_gradientn(colours=colours,limits=limits,name=legend)+scale_colour_manual(values=palette)+
      scale_x_continuous(breaks=seq_along(ordered),labels=ordered,expand=c(0,0))+
      scale_y_continuous(breaks=seq_along(ordered),labels=rev(ordered),expand=c(0,0))+
      coord_fixed()+theme_classic()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=8),
                                         axis.text.y=element_text(size=8))+
      labs(title=title,x=NULL,y=NULL,colour="Group",
           caption="Average-linkage sample-distance ordering; strip shows comparison groups.")
  }
  distance_plot <- heatmap_plot(distance,"Sample distance | all filtered VST genes","RMS distance",
                                c("#F7FBFF","#6BAED6","#08306B"))
  correlation_plot <- heatmap_plot(correlation,"Sample correlation | up to 500 variable VST genes","Pearson r",
                                   c("#B2182B","#F7F7F7","#2166AC"),c(-1,1))
  size_factors <- sizeFactors(dds)
  # tximport stores gene/sample normalization offsets, not scalar size factors.
  if(is.null(size_factors)) size_factors <- rep(NA_real_,length(samples))
  summary <- data.frame(sample_id=samples,group=groups,assigned_counts=colSums(raw),
                        genes_ge10=colSums(raw>=10),size_factor=size_factors,row.names=NULL)
  summary$normalization <- if(is.null(normalizationFactors(dds))) "Scalar size factors" else "Gene/sample offsets (tximport)"
  metrics <- rbind(data.frame(sample=samples,group=groups,metric="Gene-assigned counts (millions)",value=summary$assigned_counts/1e6),
                   data.frame(sample=samples,group=groups,metric="Genes with at least 10 counts",value=summary$genes_ge10),
                   data.frame(sample=samples,group=groups,metric="DESeq2 size factor",value=summary$size_factor))
  bar_plot <- function(frame,title,caption) {
    frame$sample <- factor(frame$sample,levels=samples)
    ggplot(frame,aes(sample,value,fill=group))+geom_col(na.rm=TRUE)+facet_wrap(~metric,scales="free_y",ncol=1)+
      scale_fill_manual(values=palette)+theme_classic()+theme(axis.text.x=element_text(angle=45,hjust=1,size=8))+
      labs(title=title,x=NULL,y=NULL,fill="Group",caption=caption)
  }
  count_plot <- bar_plot(metrics,"Count depth, detected genes and normalization",
                         "Counts are gene-assigned, not FASTQ depth. Salmon uses gene/sample offsets: scalar size factors are not applicable.")
  norm <- log2(counts(dds,normalized=TRUE)+1)
  distributions <- data.frame(sample=factor(rep(samples,each=nrow(norm)),levels=samples),
                               group=rep(groups,each=nrow(norm)),value=as.vector(norm))
  distribution_plot <- ggplot(distributions,aes(sample,value,fill=group))+geom_boxplot(outlier.shape=NA)+
    scale_fill_manual(values=palette)+theme_classic()+theme(axis.text.x=element_text(angle=45,hjust=1,size=8))+
    labs(title="Normalized expression distributions",x=NULL,y="log2(normalized count + 1)",fill="Group",
         caption="Same filtered genes in every sample; distributions are not sample-outlier tests.")
  notice <- function(title,text) ggplot()+annotate("text",x=0,y=0,label=text,size=4)+
    theme_void()+labs(title=title)
  # Missing processing QC never requests upstream files; state availability.
  summary$processing_qc_status <- "Unavailable"
  qc_plot <- notice("Retained read-processing QC","No retained processing QC supplied.\nNo reads or alignment will be regenerated.")
  files <- if(is.null(opt[["qc-files"]])) character() else fromJSON(opt[["qc-files"]])
  if(length(files)) {
    qc <- do.call(rbind,lapply(files,function(path) read.delim(gzfile(path),check.names=FALSE)))
    if(anyDuplicated(qc$run_id)) stop("Duplicate processing-QC run IDs")
    collapse <- unlist(pf$collapse)
    if(is.null(collapse)) collapse <- setNames(samples,samples)
    keys <- intersect(c("hisat2_overall_alignment","featurecounts_assigned_fraction","sense_fraction"),names(qc))
    frames <- list()
    for(i in seq_along(samples)) {
      runs <- names(collapse)[collapse == samples[i]]
      rows <- match(runs,qc$run_id)
      valid <- rows[!is.na(rows)]
      summary$processing_qc_status[i] <- if(!length(valid)) "Unavailable" else if(anyNA(rows)) "Partial" else "Available"
      for(key in keys) {
        value <- if(length(valid)) mean(qc[[key]][valid],na.rm=TRUE) else NA_real_
        frames[[length(frames)+1]] <- data.frame(sample=samples[i],group=groups[i],metric=key,value=value)
        summary[i,key] <- value
      }
      if("inferred_strandedness" %in% names(qc))
        summary[i,"inferred_strandedness"] <- paste(unique(qc$inferred_strandedness[valid]),collapse=",")
    }
    if(length(frames)) qc_plot <- bar_plot(do.call(rbind,frames),"Retained read-processing QC",
      "Means across processing runs per biological replicate; missing values are not zeros. Low sense fraction can be expected for first-strand libraries.")
  }
  summary$cooks_genes_above_cutoff <- NA_integer_
  summary$cooks_q99 <- NA_real_
  summary$cooks_status <- "Disabled"
  cooks_plot <- notice("Gene-level influence diagnostics","Cook's diagnostics disabled by configuration.")
  if(!identical(opt$cooks,"false")) {
    sizes <- table(groups)
    if(length(sizes)!=2 || any(sizes<2) || identical(pf$inference_supported,FALSE)) {
      summary$cooks_status <- "Not supported by comparison design"
      cooks_plot <- notice("Gene-level influence diagnostics","Comparison preflight or replication does not support a two-group fit.\nDescriptive plots remain available.")
    } else {
      fit <- tryCatch(DESeq(dds,fitType="mean",minReplicatesForReplace=Inf,quiet=TRUE),error=function(e)e)
      if(inherits(fit,"error")) {
        summary$cooks_status <- paste("Fit unavailable:",conditionMessage(fit))
        cooks_plot <- notice("Gene-level influence diagnostics",paste("Fit unavailable:",conditionMessage(fit)))
      } else {
        cooks <- assays(fit)[["cooks"]]
        parameters <- ncol(model.matrix(design(dds),as.data.frame(colData(dds))))
        cutoff <- qf(.99,parameters,length(samples)-parameters)
        summary$cooks_genes_above_cutoff <- colSums(cooks>cutoff,na.rm=TRUE)
        summary$cooks_q99 <- apply(cooks,2,quantile,.99,na.rm=TRUE)
        summary$cooks_status <- "Available"
        summary$cooks_reference_cutoff <- cutoff
        frame <- rbind(data.frame(sample=samples,group=groups,metric="Genes above Cook's reference cutoff",value=summary$cooks_genes_above_cutoff),
                        data.frame(sample=samples,group=groups,metric="99th percentile gene-level Cook's distance",value=summary$cooks_q99))
        cooks_plot <- bar_plot(frame,"Influence on the two-group expression model",
          sprintf("Gene-level diagnostics, not sample p-values. F(%d,%d) 99%% reference = %.2f; no count replacement.",parameters,length(samples)-parameters,cutoff))
      }
    }
  }
  pages <- list(original_pca,pca_views,distance_plot,correlation_plot,
                 arrangeGrob(count_plot,distribution_plot,ncol=1,heights=c(2,1)),qc_plot,cooks_plot)
  report <- marrangeGrob(pages,nrow=1,ncol=1,top=NULL)
  ggsave(paste0(opt$out,".pdf"),plot=report,width=12,height=9,units="in")
  write.table(summary,gzfile(paste0(opt$out,".summary.tsv.gz")),sep="\t",row.names=FALSE,quote=FALSE)
}
