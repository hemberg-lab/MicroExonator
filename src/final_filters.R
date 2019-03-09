#!/usr/bin/env Rscript

#' ---
#' title: Micro-exon final filtering report
#' author: Guillermo Parada
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
#' 
#' ## Loading libraries
#' 
#' The libraries which are used by this script are:"
#' *ggplot2
#' *reshape2
#' *stringr
#' *mixtools
#' *simecol
#' *data.table
#' *optparse


library(ggplot2)
library(reshape2)
library(stringr)
library(mixtools)
library(simecol)
library(data.table)
#library(optparse)


#' ## Input micro-exon profiling
#' 
#' De-novo discovery of micro-exons by uExonator relies on the detection of inserted sequenses over
#' exon-exon junctions, which then are re-mapped inside the cognate introns. Inserted sequences smaller
#' than 6 nucleotides are very likely to be mapped by chance, therefore detected micro-exons smaller than
#' 6 nt are prone to be artefacts by sequencing error or genomic variations. The following plot shows the
#' spurious micro-exon/intron match probability distribution for micro-exon in between 1-15 nt.




#option_list <- list(
#  make_option(c("-met", "--micro_exon_table"), type="character", default=NULL, 
#              help="Micro-exon centric table", metavar="character"),
#  make_option(c("-c", "--micro_exon_coverages"), type="character", default=NULL, 
#              help="Micro-exon coverage table", metavar="character"),  
#  make_option(c("-o", "--out"), type="character", default="out.txt", 
#              help="output file name [default= %default]", metavar="character")
#); 

#opt_parser <- OptionParser(option_list=option_list);
#opt <- parse_args(opt_parser);



ME_centric_raw <- read.delim("~/Google_Drive/Results/ME/Single_cell/TOTAL.sam.row_ME.filter1.ME_centric", header=FALSE, stringsAsFactors=FALSE)
colnames(ME_centric_raw) <- c('ME', 'transcript', 'sum_total_coverage', 'total_SJs', 'total_coverages', 'len_micro_exon_seq_found', 'micro_exon_seq_found', 'total_number_of_micro_exons_matches', 'U2_scores', 'mean_conservations_vertebrates', 'P_MEs', 'total_ME')

ME_centric_raw <- data.table(ME_centric_raw)

ggplot(ME_centric_raw[len_micro_exon_seq_found<=15, ],
       aes(x=factor(len_micro_exon_seq_found), y=P_MEs) ) +
  geom_violin(scale = "width") +
  xlab("Micro-exon leght") +
  ylab("Spurious micro-exon/intron match probability") +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))

#' The higher spurious micro-exon/intron match probability is reflected on the number of micro-exon/intron
#' matches inside

ggplot(ME_centric_raw[len_micro_exon_seq_found<=15, ],
       aes(x=factor(len_micro_exon_seq_found), y=total_number_of_micro_exons_matches) ) + 
  geom_jitter() +
  ylim(0,100) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))


#'     True 
#' splicing events relies on splicing signals, therefore false micro-exons will have weaker splicing signals
#' than the true micro-exons. The following plot show the distribution of U2/GT-AG splicing signal strengh
#' (U2_score) for population the total micro-exons and longer or equal then 3, 6, and 9 nt. Micro-exons equal
#' or longer than 9 nt are less prone to be artefacts, therefore have a U2_score distribution which is expected
#' from real splicng events.


#ME_matches <- unlist(strsplit(ME_centric_raw$total_ME, "[,]"))
#ME_matches <- read.table(text=ME_matches, sep="|")
#colnames(ME_matches) <- c("ME", "U2_score", "Vertebrate_conservation", "Primate_conservation")
#ME_matches$Filter = "None"
#ME_centric_raw_longer_3 <- subset(ME_centric_raw, len_micro_exon_seq_found>=3)
#ME_centric_raw_longer_6 <- subset(ME_centric_raw, len_micro_exon_seq_found>=6)
#ME_centric_raw_longer_9 <- subset(ME_centric_raw, len_micro_exon_seq_found>=9)
#ME_matches_3 <- unlist(strsplit(ME_centric_raw_longer_3$total_ME, "[,]"))
#ME_matches_3 <- read.table(text=ME_matches_3, sep="|")
#colnames(ME_matches_3) <- c("ME", "U2_score", "Vertebrate_conservation", "Primate_conservation")
#ME_matches_3$Filter = ">=3"
#ME_matches_6 <- unlist(strsplit(ME_centric_raw_longer_6$total_ME, "[,]"))
#ME_matches_6 <- read.table(text=ME_matches_6, sep="|")
#colnames(ME_matches_6) <- c("ME", "U2_score", "Vertebrate_conservation", "Primate_conservation")
#ME_matches_6$Filter = ">=6"
#ME_matches_9 <- unlist(strsplit(ME_centric_raw_longer_9$total_ME, "[,]"))
#ME_matches_9 <- read.table(text=ME_matches_9, sep="|")
#colnames(ME_matches_9) <- c("ME", "U2_score", "Vertebrate_conservation", "Primate_conservation")
#ME_matches_9$Filter = ">=9"
#ME_matches_Filters <- rbind(ME_matches, ME_matches_3, ME_matches_6, ME_matches_9)
#ggplot(ME_matches_Filters, aes(x=U2_score, ..density.., colour=Filter)) +
#  geom_freqpoly(binwidth=5) +
#  xlim(40, 100) +
#  theme(panel.background = element_rect(fill = 'white', colour = 'black'))


R -e "rmarkdown::render('../../../Software/Micro-Exonator/src/final_filters.R',output_file='output.html')"