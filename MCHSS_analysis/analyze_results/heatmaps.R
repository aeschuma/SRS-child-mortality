# Austin Schumacher
# 2/13/2020
# Make heatmaps of the change in mortality (1996 to 2015) from final estimates from model fit to MCHSS data

rm(list=ls())

## load libraries
library(INLA); library(scales); library(RColorBrewer); library(tidyverse); library(ggpubr);
library(grid); library(gridExtra); library(viridis);

#---------------
## define directories
#---------------

# directory for results
savedir <- "../../../../Dropbox/SRS-child-mortality-output/"

# create folders to store results if necessary
if (!file.exists(paste0(savedir, "graphs"))) {
    dir.create(paste0(savedir, "graphs"))
}
if (!file.exists(paste0(savedir, "graphs/heatmaps"))) {
    dir.create(paste0(savedir, "graphs/heatmaps"))
}

# load our results
mod_inla <- readRDS(paste0(savedir, "results/china_results_pcpriors.RDS"))
results <- readRDS(paste0(savedir, "results/china_results_wide.RDS"))
res <- readRDS(paste0(savedir, "results/china_results_long.RDS"))

#-----------------
## Various posterior predictions to calculate posterior median change in mortality
#-----------------

# data dimensions
nyear <- length(unique(results$year))
ncause <- length(unique(results$cause))
nreg <- length(unique(results$reg2))
nageg <- length(unique(results$agegp))
years <- unique(results$year)

# number of samples to take for posterior sampling
nsamples <- 1000

# predictions for log mortaity rate
samples <- inla.posterior.sample(nsamples,mod_inla)

contents <- mod_inla$misc$configs$contents

id.pred <- which(contents$tag=="Predictor")
ind.pred <- contents$start[id.pred]-1 + (1:contents$length[id.pred])

samples.pred <- lapply(samples, function(x) {
  x$latent[ind.pred] - results$logpy
})
s.pred <- t(matrix(unlist(samples.pred), byrow = T, nrow = length(samples.pred)))
samples.pred <- NULL

# FE only predictions
id.fe <- which(contents$tag=="(Intercept)")
ind.fe <- contents$start[id.fe]-1 + (1:sum(contents$length[id.fe:length(contents$length)]))

X <- model.matrix(~factor(agegp)*factor(reg2) +
                    factor(agegp)*factor(cause) +
                    factor(reg2)*factor(cause), 
                  data = results)
samples.fe <- lapply(samples, function(x) {
  as.vector(X %*% x$latent[ind.fe])
})
s.fe <- t(matrix(unlist(samples.fe), byrow = T, nrow = length(samples.fe)))
samples.fe <- NULL

# RW predictions
id.rw <- which(contents$tag=="year")
ind.rw <- contents$start[id.rw]-1 + (1:contents$length[id.rw])

samples.rw <- lapply(samples, function(x) {
  x$latent[ind.rw]
})
s.rw <- t(matrix(unlist(samples.rw), byrow = T, nrow = length(samples.rw)))
samples.rw <- NULL

# RE predictions
id.re <- which(contents$tag=="factor(obs)")
ind.re <- contents$start[id.re]-1 + (1:contents$length[id.re])

samples.re <- lapply(samples, function(x) {
  x$latent[ind.re]
})
s.re <- t(matrix(unlist(samples.re), byrow = T, nrow = length(samples.re)))
samples.re <- NULL

# FE + RW only predictions
s.fe.rw <- matrix(NA,nrow=nrow(results),ncol=nsamples)
for (i in 1:nsamples) {
  if (i == 1) cat(paste("Starting sim 1 ...")); flush.console()
  if (i %% 10 == 0) cat(paste(i, "...")); flush.console()
  rwcount <- 0
  for (j in 1:length(unique(results$rw_index))) {
    for (k in 1:nyear) {
      rwcount <- rwcount + 1
      s.fe.rw[results$rw_index==j & results$year==years[k],i] <- 
        s.fe[results$rw_index==j & results$year==years[k],i] + 
        s.rw[rwcount,i]
    }
  }
}

#------------
## Heatmap plot of change in mortality for each age-region-cause strata
#------------

# make heatmap data frame to store results
heatmap <- results[!duplicated(results$agereg2cause),c("agegp_name","reg2","cause_name","agereg2cause")]
heatmap <- heatmap[order(heatmap$agegp_name,heatmap$reg2,heatmap$cause_name),]
heatmap$diff_fe_rw_m <- NA
heatmap$diff_fe_rw_l <- NA
heatmap$diff_fe_rw_u <- NA
heatmap$rel_diff_fe_rw_m <- NA
heatmap$rel_diff_fe_rw_l <- NA
heatmap$rel_diff_fe_rw_u <- NA
for (i in 1:length(unique(results$agereg2cause))) {
  draws_tmp2 <- (exp(s.fe.rw[results$agereg2cause == heatmap$agereg2cause[i] & results$year_name == 2015,]) - 
                   exp(s.fe.rw[results$agereg2cause == heatmap$agereg2cause[i] & results$year_name == 1996,]))
  heatmap$diff_fe_rw_m[i] <- median(draws_tmp2)
  heatmap$diff_fe_rw_l[i] <- quantile(draws_tmp2,0.1)
  heatmap$diff_fe_rw_u[i] <- quantile(draws_tmp2,0.9)
  
  draws_tmp4 <- (exp(s.fe.rw[results$agereg2cause == heatmap$agereg2cause[i] & results$year_name == 2015,]) - 
                   exp(s.fe.rw[results$agereg2cause == heatmap$agereg2cause[i] & results$year_name == 1996,])) / 
    exp(s.fe.rw[results$agereg2cause == heatmap$agereg2cause[i] & results$year_name == 1996,])
  heatmap$rel_diff_fe_rw_m[i] <- median(draws_tmp4)*100
  heatmap$rel_diff_fe_rw_l[i] <- quantile(draws_tmp4,0.1)*100
  heatmap$rel_diff_fe_rw_u[i] <- quantile(draws_tmp4,0.9)*100
}

## set names
cause_names <- unique(heatmap$cause_name)
agegp_names <- unique(heatmap$agegp_name)
reg2_names <- unique(heatmap$reg2)

## make heatmap long data frame
heatmap_long <- gather(heatmap, measure, value, diff_fe_rw_m:rel_diff_fe_rw_u)

# settings
textsize <- 22

ggplot(data = heatmap_long[heatmap_long$measure == "diff_fe_rw_m" 
                           & heatmap_long$agegp_name == "0-6d",],
       aes(x = cause_name, y = reg2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(high = "white",
                      low = muted("green"),
                      breaks=c(-0.8,-0.4,0),
                      limits = c(-0.8,0)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle(NULL) +
  theme_light() +
  theme(legend.position = "top",
        legend.key.size = unit(1.8, "cm"),
        legend.text = element_text(size = textsize, angle = 45, vjust = 0.5, face = "bold"),
        legend.title = element_blank(),
        plot.title = element_blank(),
        strip.text = element_text(size = textsize, face = "bold"),
        axis.text.y = element_text(size = textsize, face = "bold"),
        axis.text.x = element_text(size = textsize, angle = 45, face = "bold", vjust = 1, hjust = 1))
ggsave(filename = paste0(savedir, "graphs/heatmaps/mx_pcpriors_fe_rw_m_heatmap_0-6d.pdf"), 
       device = "pdf", width = 7, height = 8)

ggplot(data = heatmap_long[heatmap_long$measure == "rel_diff_fe_rw_m" 
                           & heatmap_long$agegp_name == "1-5m",],
       aes(x = cause_name, y = reg2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(high = muted("red"), mid = "white",
                      low = muted("blue"),
                      breaks=c(-100,-50,0,50,100),
                      limits = c(-100,100)) +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle(NULL) +
  theme_light() +
  theme(legend.position = "top",
        legend.key.size = unit(1.8, "cm"),
        legend.text = element_text(size = textsize, angle = 45, vjust = 0.5, face = "bold"),
        legend.title = element_blank(),
        plot.title = element_blank(),
        strip.text = element_text(size = textsize, face = "bold"),
        axis.text.y = element_text(size = textsize, face = "bold"),
        axis.text.x = element_text(size = textsize, angle = 45, face = "bold", vjust = 1, hjust = 1))
ggsave(filename = paste0(savedir, "graphs/heatmaps/mx_pcpriors_rel_fe_rw_m_heatmap_1-5m.pdf"), 
       device = "pdf", width = 7, height = 8)
