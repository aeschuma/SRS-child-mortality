## Austin Schumacher
## 1/3/2019
## Simulation to motivate our method 
## Goals
## - demonstrate problems with two-stage estimation framework
## - explore how a simple two-stage framework breaks in certain situations
## - compare this with a unified framework
## This code:
## - compiles results from 'simulation_comparison_no_correlation.R'
## - plots comparisons of the multistage and unified models

rm(list=ls())

## load libraries
library(scales); library(RColorBrewer); library(ggplot2); library(cowplot);

# directory for results
savedir <- "../../../Dropbox/SRS-child-mortality-output/"

# create folders to store results if necessary
# note: the results folders were created in the simulation code already
if (!file.exists(paste0(savedir, "graphs"))) {
    dir.create(paste0(savedir, "graphs"))
}
if (!file.exists(paste0(savedir, "graphs/sims"))) {
    dir.create(paste0(savedir, "graphs/sims"))
}
if (!file.exists(paste0(savedir, "graphs/sims/no_correlation"))) {
    dir.create(paste0(savedir, "graphs/sims/no_correlation"))
}

## parameters
multinomial <- FALSE
nsim <- 100
ncause <- 2
nstrata <- 6*6*20
exposures <- c(1000,10000,100000)
nexposures <- length(exposures)

## store results
results_2stage <- array(data=NA,dim=c(nexposures,nsim,4))
results_casm <- array(data=NA,dim=c(nexposures,nsim,4))

## load files
fl_casm <- grep(paste0(savedir, "results/sims/tmp/no_correlation/motive_res_casm_iid_intonly_casm_multinom_",multinomial,"_ncause_",ncause),list.files(),value=T)
fl_2stage <- grep(paste0(savedir, "results/sims/tmp/no_correlation/motive_res_iid_intonly_multinom_",multinomial,"_ncause_",ncause),list.files(),value=T)

for (i in 1:length(fl_2stage)) {

    sim_2stage <- as.numeric(substr(fl_2stage[i],
                                    gregexpr(pattern ='_',fl_2stage[i])[[1]][9]+1,
                                    gregexpr(pattern ='_',fl_2stage[i])[[1]][10]-1))
    exposure_2stage <- as.numeric(substr(fl_2stage[i],
                                         gregexpr(pattern ='_',fl_2stage[i])[[1]][11]+1,
                                         gregexpr(pattern ='.rds',fl_2stage[i])[[1]][1]-1))

    results_2stage[which(exposures==exposure_2stage),sim_2stage,] <- readRDS(fl_2stage[i])
    
    sim_casm <- as.numeric(substr(fl_casm[i],
                                  gregexpr(pattern ='_',fl_casm[i])[[1]][11]+1,
                                  gregexpr(pattern ='_',fl_casm[i])[[1]][12]-1))
   
    exposure_casm <- as.numeric(substr(fl_casm[i],
                                       gregexpr(pattern ='_',fl_casm[i])[[1]][13]+1,
                                       gregexpr(pattern ='.rds',fl_casm[i])[[1]][1]-1))
    
    results_casm[which(exposures==exposure_casm),sim_casm,] <- readRDS(fl_casm[i])
}

## save results
saveRDS(results_2stage,file=paste0(savedir, "results/sims/res_motive_iid_intonly_multinom_",multinomial,"_ncause_",ncause,"_2stage_",Sys.Date(),
                                   "_nsim_",nsim,"_max_exp_",max(exposures),"_nstrata_",nstrata,".rds"))
saveRDS(results_casm,file=paste0(savedir, "results/sims/res_motive_iid_intonly_multinom_",multinomial,"_ncause_",ncause,"_rw_casm_",Sys.Date(),
                                 "_nsim_",nsim,"_max_exp_",max(exposures),"_nstrata_",nstrata,".rds"))

# delete all temporary sim data (to save space)
if (file.exists(paste0(savedir, "results/sims/res_motive_iid_intonly_multinom_",multinomial,"_ncause_",ncause,"_2stage_",Sys.Date(),
                       "_nsim_",nsim,"_max_exp_",max(exposures),"_nstrata_",nstrata,".rds"))) {
    unlink(unlist(fl_2stage))
}
if (file.exists(paste0(savedir, "results/sims/res_motive_iid_intonly_multinom_",multinomial,"_ncause_",ncause,"_rw_casm_",Sys.Date(),
                       "_nsim_",nsim,"_max_exp_",max(exposures),"_nstrata_",nstrata,".rds"))) {
    unlink(unlist(fl_casm))
}

# load results (this is just here for convenience---only necessary if we're redoing the plots)
# results_2stage <- readRDS(paste0(savedir, "results/sims/res_motive_iid_intonly_multinom_FALSE_ncause_2_2stage_2020-01-28_ncause_2_nsim_100_max_exp_1e+05_nstrata_720.rds"))
# results_casm <- readRDS(paste0(savedir, "results/sims/res_motive_iid_intonly_multinom_FALSE_ncause_2_rw_casm_2020-01-28_ncause_2_nsim_100_max_exp_1e+05_nstrata_720.rds"))

# aggregate 
res_2stage_agg <- apply(results_2stage,c(1,3),mean,na.rm=T)
res_casm_agg <- apply(results_casm,c(1,3),mean,na.rm=T)

# create data frame
plot.df <- data.frame(model = rep(c(rep("multistage", nrow(res_2stage_agg)), 
                                rep("unified", nrow(res_casm_agg))), 4),
                      exposure = rep(rep(exposures,2), 4),
                      measure = rep(c("Relative bias","Absolute bias","Coverage","CI width"),
                                    each = length(exposures)*2),
                      value = c(c(res_2stage_agg[,1], res_casm_agg[,1]), 
                                c(res_2stage_agg[,2], res_casm_agg[,2]),
                                c(res_2stage_agg[,3], res_casm_agg[,3]),
                                c(res_2stage_agg[,4], res_casm_agg[,4])))
plot.df$model <- as.character(plot.df$model)

###############
# PLOTS:
#   - for CSMR
#       - bias, rel bias, coverage, CI width
###############

plot.list <- list()
i <- 0
theme_set(theme_cowplot())
for (outcome in unique(plot.df$measure)) {
    i <- i + 1
    if (outcome %in% c("Relative bias","Absolute bias","CI width")) ref_yint <- 0
    if (outcome == "Coverage") ref_yint <- 0.95
    p <- ggplot(data = plot.df[plot.df$measure == outcome,], 
                aes(x = exposure/1000, y = value,
                    col = model, linetype = model)) + 
        geom_line(size = 1.25) +
        geom_hline(yintercept = ref_yint, alpha = 0.4, size = 0.75) + 
        xlab('Exposure (1000s)') +
        ylab(outcome) +
        theme(axis.text = element_text(size=14),
              axis.title = element_text(size=16,face="bold"),
              legend.text = element_text(size=13),
              legend.title = element_text(size=16,face="bold"))
    plot.list[[i]] <- p
}

prow <- plot_grid(plot.list[[1]] + theme(legend.position="none"), 
                  plot.list[[3]] + theme(legend.position="none"), 
                  plot.list[[4]] + theme(legend.position="none"),
                  labels = "AUTO", nrow = 1)
legend <- get_legend(plot.list[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(prow, legend, rel_widths = c(3, .4))
ggsave(paste0(savedir, "graphs/sims/no_correlation/res_motive_iid_intonly_multinomial_",multinomial,"_ncause_",ncause,"_nsim_",nsim,"_max_exposure_",max(exposures),"_nstrata_",nstrata,".pdf"),
       width=11,height=3)

## plot with no bias
prow2 <- plot_grid(plot.list[[3]] + theme(legend.position="none"), 
                   plot.list[[4]] + theme(legend.position="none"),
                   labels = "AUTO", nrow = 1)
legend2 <- get_legend(plot.list[[3]] + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(prow2, legend2, rel_widths = c(2, .3))
ggsave(paste0(savedir, "graphs/sims/no_correlation/res_motive_nobias_iid_intonly_multinom_",multinomial,"_ncause_",ncause,"_nsim_",nsim,"_max_exposure_",max(exposures),"_nstrata_",nstrata,".pdf"),
       width=12,height=3)
