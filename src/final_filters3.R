log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library(ggplot2)
library(reshape2)
library(stringr)
library(mixtools)
#library(simecol)
library(data.table)



contrast <- c("condition", snakemake@params[["contrast"]])



ME_centric_raw <- data.table(ME_centric_raw)

ME_number_files_detected <- fread(ME_count_round2)

ME_centric_raw_filter_uniq <- ME_centric_raw[ME %in% ME_number_files_detected[N_samples >=min_number_files_detected, ME], ]
uniq_seq_filter <-  ME_centric_raw_filter_uniq[, ME]


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
