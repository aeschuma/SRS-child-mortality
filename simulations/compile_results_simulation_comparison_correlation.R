## Austin Schumacher
## 1/17/2019
## Simulation to motivate our method 
## Goals
## - demonstrate problems with two-stage estimation framework
## - explore how a simple two-stage framework breaks in certain situations
## - compare this with a unified framework
## This code:
## - compiles results from 'simulation_comparison_correlation.R'
## - plots comparisons of the multistage and unified models

rm(list=ls())

## set the root depending on operating system
root <- ifelse(Sys.info()[1]=="Darwin","~/",
               ifelse(Sys.info()[1]=="Windows","P:/",
                      ifelse(Sys.info()[1]=="Linux","/home/students/aeschuma/",
                             stop("Unknown operating system"))))

## the following code makes rstan work on the Box server
if (root == "P:/") {
    Sys.setenv(HOME="C:/Users/aeschuma",
               R_USER="C:/Users/aeschuma",
               R_LIBS_USER="C:/Users/aeschuma/R_libraries")
    .libPaths("C:/Users/aeschuma/R_libraries")
} else if (root == "/home/students/aeschuma/") {
    Sys.setenv(HOME=root,
               R_USER=root,
               R_LIBS_USER=paste0(root,"R/x86_64-pc-linux-gnu-library/3.5"))
    .libPaths(paste0(root,"R/x86_64-pc-linux-gnu-library/3.5"))
}

## load libraries
library(scales); library(RColorBrewer); library(ggplot2); library(cowplot);

## define directories

# working directory for code
wd <- paste(root,"Desktop/dissertation/motivation_sim",sep="")

# directory to save results
savedir <- paste(root,"Dropbox/dissertation_2/cause_specific_child_mort/motivation_sim",sep="")

## parameters
intercept_only <- FALSE
linear_time <- FALSE
time_rw <- FALSE
interactions <- FALSE
multinomial <- FALSE
nsim <- 100
ncause <- 2
rhos <- c(-0.5,0,0.5)
sigma_res <- c(0.01,0.1,1)

## store results
results_2stage <- array(data = NA, dim = c(length(sigma_res), length(rhos),nsim,4), 
                        dimnames = list(sigma = sigma_res,
                                        rho = rhos,
                                        sim = 1:nsim,
                                        measure = c("Relative bias","Absolute bias","Coverage","CI width")))
results_casm <- array(data = NA, dim = c(length(sigma_res),length(rhos),nsim,4),
                      dimnames = list(sigma = sigma_res,
                                      rho = rhos,
                                      sim = 1:nsim,
                                      measure = c("Relative bias","Absolute bias","Coverage","CI width")))

## set directory
setwd(paste(savedir,"results",sep="/"))

## load files
fl_2stage <- grep(paste0("motive_corr_res",
                       "_intonly_",intercept_only,
                       "_lineartime_",linear_time,
                       "_timerw_",time_rw,
                       "_interactions_",interactions,
                       "_multinom_",multinomial,
                       "_ncause_",ncause),
                  list.files(),
                  value=T)
fl_casm <- grep(paste0("motive_corr_res_casm",
                       "_intonly_",intercept_only,
                       "_lineartime_",linear_time,
                       "_timerw_",time_rw,
                       "_interactions_",interactions,
                       "_multinom_",multinomial,
                       "_ncause_",ncause),
                list.files(),
                value=T)

for (i in 1:length(fl_2stage)) {
    
    # example file name
    # motive_corr_res_casm_intonly_FALSE_lineartime_TRUE_timerw_FALSE_interactions_FALSE_multinom_FALSE_ncause_2_sim_88_sigma_0.01_rho_-0.5.rds"
    
    sim_2stage <- as.numeric(substr(fl_2stage[i],
                                    gregexpr(pattern ='_',fl_2stage[i])[[1]][16]+1,
                                    gregexpr(pattern ='_',fl_2stage[i])[[1]][17]-1))
    sigma_2stage <- as.numeric(substr(fl_2stage[i],
                                         gregexpr(pattern ='_',fl_2stage[i])[[1]][18]+1,
                                         gregexpr(pattern ='_',fl_2stage[i])[[1]][19]-1))
    rho_2stage <- as.numeric(substr(fl_2stage[i],
                                       gregexpr(pattern ='_',fl_2stage[i])[[1]][20]+1,
                                       gregexpr(pattern ='.rds',fl_2stage[i])[[1]][1]-1))
    
    results_2stage[as.numeric(which(sigma_res==sigma_2stage)),
                   as.numeric(which(rhos==rho_2stage)),
                   sim_2stage,] <- readRDS(fl_2stage[i])
    
    sim_casm <- as.numeric(substr(fl_casm[i],
                                    gregexpr(pattern ='_',fl_casm[i])[[1]][17]+1,
                                    gregexpr(pattern ='_',fl_casm[i])[[1]][18]-1))
    sigma_casm <- as.numeric(substr(fl_casm[i],
                                      gregexpr(pattern ='_',fl_casm[i])[[1]][19]+1,
                                      gregexpr(pattern ='_',fl_casm[i])[[1]][20]-1))
    rho_casm <- as.numeric(substr(fl_casm[i],
                                  gregexpr(pattern ='_',fl_casm[i])[[1]][21]+1,
                                  gregexpr(pattern ='.rds',fl_casm[i])[[1]][1]-1))
    
    results_casm[as.numeric(which(sigma_res==sigma_casm)),
                 as.numeric(which(rhos==rho_casm)),
                 sim_casm,] <- readRDS(fl_casm[i])
}

## save results
saveRDS(results_2stage,file=paste0("res_corr_motive_2stage_",Sys.Date(),
                                   "_intonly_",intercept_only,
                                   "_lineartime_",linear_time,
                                   "_timerw_",time_rw,
                                   "_interactions_",interactions,
                                   "_multinom_",multinomial,
                                   "_ncause_",ncause,
                                   "_nsim_",nsim,
                                   ".rds"))
saveRDS(results_casm,file=paste0("res_corr_motive_casm_",Sys.Date(),
                                 "_intonly_",intercept_only,
                                 "_lineartime_",linear_time,
                                 "_timerw_",time_rw,
                                 "_interactions_",interactions,
                                 "_multinom_",multinomial,
                                 "_ncause_",ncause,
                                 "_nsim_",nsim,
                                 ".rds"))

## delete all tmp results
if (file.exists(paste0("res_corr_motive_2stage_",Sys.Date(),
                       "_intonly_",intercept_only,
                       "_lineartime_",linear_time,
                       "_timerw_",time_rw,
                       "_interactions_",interactions,
                       "_multinom_",multinomial,
                       "_ncause_",ncause,
                       "_nsim_",nsim,
                       ".rds"))) {
    unlink(unlist(fl_2stage))
}
if (file.exists(paste0("res_corr_motive_casm_",Sys.Date(),
                       "_intonly_",intercept_only,
                       "_lineartime_",linear_time,
                       "_timerw_",time_rw,
                       "_interactions_",interactions,
                       "_multinom_",multinomial,
                       "_ncause_",ncause,
                       "_nsim_",nsim,
                       ".rds"))) {
    unlink(unlist(fl_casm))
}

# load results (only if redoing plots)
results_2stage <- readRDS("res_corr_motive_2stage_2020-01-20_intonly_FALSE_lineartime_FALSE_timerw_FALSE_interactions_FALSE_multinom_FALSE_ncause_2_nsim_100.rds")
results_casm <- readRDS("res_corr_motive_casm_2020-01-20_intonly_FALSE_lineartime_FALSE_timerw_FALSE_interactions_FALSE_multinom_FALSE_ncause_2_nsim_100.rds")

# aggregate 
res_2stage_agg <- apply(results_2stage,c(1,2,4),mean,na.rm=T)
res_casm_agg <- apply(results_casm,c(1,2,4),mean,na.rm=T)

# reshape and combine
plot_2stage_tbl <- res_2stage_agg %>%
    as.tbl_cube(met_name = "value") %>%
    as_tibble
plot_2stage_tbl$model <- "multistage"
plot_casm_tbl <- res_casm_agg %>%
    as.tbl_cube(met_name = "value") %>%
    as_tibble
plot_casm_tbl$model <- "unified"
plot.df <- as.data.frame(rbind(plot_2stage_tbl,plot_casm_tbl))

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
                aes(x = sigma, y = value, linetype = factor(rho),
                    col = model, pch = model)) + 
        geom_line(size = 0.75) +
        geom_point(size = 3) + 
        geom_hline(yintercept = ref_yint, alpha = 0.4, size = 0.75) + 
        xlab(expression(sigma^2)) +
        ylab(outcome) + 
        labs(linetype = expression(rho)) +
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
ggsave(paste0("../graphs/res_corr_motive_",Sys.Date(),
              "_intonly_",intercept_only,
              "_lineartime_",linear_time,
              "_timerw_",time_rw,
              "_interactions_",interactions,
              "_multinom_",multinomial,
              "_ncause_",ncause,
              "_nsim_",nsim,
              ".pdf"),
       width=11,height=3)

## plot with no bias
prow2 <- plot_grid(plot.list[[3]] + theme(legend.position="none"), 
                   plot.list[[4]] + theme(legend.position="none"),
                   labels = "AUTO", nrow = 1)
legend2 <- get_legend(plot.list[[3]] + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(prow2, legend2, rel_widths = c(2, .3))
ggsave(paste0("../graphs/res_corr_motive_nobias_",Sys.Date(),
              "_intonly_",intercept_only,
              "_lineartime_",linear_time,
              "_timerw_",time_rw,
              "_interactions_",interactions,
              "_multinom_",multinomial,
              "_ncause_",ncause,
              "_nsim_",nsim,
              ".pdf"),
       width=12,height=3)
