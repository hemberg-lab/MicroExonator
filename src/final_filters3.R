log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(ggplot2)
library(reshape2)
library(stringr)
library(mixtools)
#library(simecol)
library(data.table)






plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}
plot_mix_comps_total <- function(x, mu1, sigma1, lam1, mu2, sigma2, lam2) {
 lam1* dnorm(x, mu1, sigma1) + lam2* dnorm(x, mu2, sigma2)
}
ggplot_mix_comps <-function(mixmdl, title) {
data.frame(x = mixmdl$x) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black",
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                colour = "green", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps_total,
                args = list(mixmdl$mu[1], mixmdl$sigma[1], mixmdl$lambda[1],  mixmdl$mu[2], mixmdl$sigma[2], mixmdl$lambda[2] ),
                colour = "black", lwd = 0.5, lty=2) +
     ggtitle(title)+
  ylab("Density") +
  xlab("U2 Score") +
   theme(panel.background = element_rect(fill = 'white', colour = 'black'))
}




ME_centric_raw <- fread(snakemake@params[["ME_table"]] )
colnames(ME_centric_raw) <- c('ME', 'transcript', 'sum_total_coverage', 'total_SJs', 'total_coverages', 'len_micro_exon_seq_found', 'micro_exon_seq_found', 'total_number_of_micro_exons_matches', 'U2_scores', 'mean_conservations_vertebrates', 'P_MEs', 'total_ME')

ME_number_files_detected <- fread(snakemake@params[["ME_count_round2"]])




ME_centric_raw_filter_uniq <- ME_centric_raw[ME %in% ME_number_files_detected[N_samples >=min_number_files_detected, ME], ]
uniq_seq_filter <-  ME_centric_raw_filter_uniq[, ME]


#ME_number_files_detected <- ME_count_round2[sum_ME_SJ_coverage_up_down_uniq>=5, .N, by=ME]
ME_centric_raw_filtered <- ME_centric_raw[ME %in% ME_number_files_detected[N_samples >=min_number_files_detected, ME] &  len_micro_exon_seq_found>=3,]
uniq_seq_filter <-  ME_centric_raw_filtered[, ME]
ME_matches_filter <- ME_matches[ME %in% uniq_seq_filter , ]
ME_matches_filter <- ME_matches_filter[sample(dim(ME_matches_filter)[1])]
ME_matches_filter <- unique(ME_matches_filter, by = "ME_max_U2")
fit_U2_score <- normalmixEM(ME_matches_filter$U2_score, maxit = 10000, epsilon = 1e-05)
ggplot_mix_comps(fit_U2_score, "Mixture model Micro-exon >=3 after coverge filter")
post.df <- as.data.frame(cbind(x = fit_U2_score$x, fit_U2_score$posterior))




ME_final <- ME_centric_raw[ME %in% ME_number_files_detected[N >=min_number_files_detected, ME] & len_micro_exon_seq_found>=3, ]
if(fit_U2_score$mu[1]<=fit_U2_score$mu[2]){
  ME_final$ME_P_value <-  1 - (1 - approx(post.df$x, post.df$comp.1, ME_final$U2_scores)$y * ME_final$P_MEs) / ME_final$total_number_of_micro_exons_matches
} else {
  ME_final$ME_P_value <-  1 - (1 - approx(post.df$x, post.df$comp.2, ME_final$U2_scores)$y * ME_final$P_MEs)/ ME_final$total_number_of_micro_exons_matches
}



error_sim_asim <- NULL
for (i in seq(0.1,1, 0.01)){
  len_freq <- ME_final[ME_P_value <=i, .N, by=len_micro_exon_seq_found]
  sim <-  len_freq[len_micro_exon_seq_found%%3==0, ]
  asim <- len_freq[len_micro_exon_seq_found%%3!=0, ]
  sim_fit <- lm(data = sim,  N ~ len_micro_exon_seq_found)
  sim_error <- summary(sim_fit)$sigma
  asim_fit <- lm(data = asim,  N ~ len_micro_exon_seq_found)
  asim_error <- summary(asim_fit)$sigma
  error_sim_asim <- rbind(error_sim_asim, c(i, sim_error, asim_error))
}
error_sim_asim <- data.table(error_sim_asim)
colnames(error_sim_asim) <- c("P_ME", "sim_error", "asim_error")
error_sim_asim_melt <- melt(error_sim_asim, id.vars = "P_ME", value.name = "error")
ggplot(error_sim_asim, aes(x=P_ME, y=sim_error+asim_error)) +
  geom_line() +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))
error_sim_asim[ , total_error:=(sim_error + asim_error)]
P_ME_fit <- error_sim_asim[ total_error==error_sim_asim[, min(total_error)], P_ME]


ME_final[ME_P_value <= P_ME_fit, ME_type:="IN"]
ME_final[ME_P_value > P_ME_fit, ME_type:="OUT"]
ME_final[ME_P_value > P_ME_fit & mean_conservations_vertebrates>=2, ME_type:="RESCUED"]
ME_final[,ME_type:=factor(ME_type, levels=c("OUT", "RESCUED", "IN"))]
