## Austin Schumacher
## 2/29/2020
## Fit the final model in INLA on MCHSS data, and plot results

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
}

## load libraries
library(INLA); library(scales); library(RColorBrewer); library(tidyverse); library(ggpubr);
library(grid); library(gridExtra); library(viridis);

########################
## code running options
########################

#####################
## modeling options
#####################
quantiles <- c(0.025,0.1,0.5,0.9,0.975)
nsamples <- 1000

######################
## define directories
######################

# working directory for code
wd <- paste(root,"Desktop/dissertation/estimation_china",sep="")

# directory to save results
savedir <- paste(root,"Dropbox/dissertation_2/cause_specific_child_mort/estimation_china",sep="")

## set directory
setwd(savedir)

################
## read and format data
################

# load MCHSS data
chn_all <- readRDS("../china_data/chn_srs_formatted.RDS")

## format data
chn_all$logpy <- log(chn_all$exposure)
chn_all$ageRW <- ifelse(chn_all$agegp==1,1,2)
chn_all$causeRW <- as.character(chn_all$cause_name)
chn_all$causeRW[chn_all$cause_name=="diarrhea" | chn_all$cause_name=="other grp 1"] <- "diarrhea and other grp 1"
chn_all$causeRW[chn_all$cause_name=="congenital anomalies" | chn_all$cause_name=="other non-communicable"] <- "cong. anomalies and other n.c."

## cause dimensions
ncause <- length(unique(chn_all$cause))
ncause_young <- max(chn_all$cause_young_ind)
ncause_old <- max(chn_all$cause_old_ind)
young_causes <- unique(chn_all$cause_ind[chn_all$cause_young_ind != 0])
old_causes <- unique(chn_all$cause_ind[chn_all$cause_old_ind !=0 ])
nyear <- length(unique(chn_all$year))
years <- unique(chn_all$year)
nageg <- length(unique(chn_all$agegp))
young_ages <- c(1,2,3)
old_ages <- c(4,5,6)
chn_all$reg2 <- interaction(chn_all$reg,chn_all$res)

## delete the rows corresponding to old ages with causes that we think only affects young ages
chn_all <- chn_all[!(chn_all$agegp %in% old_ages & !(chn_all$cause_ind %in% old_causes)),]

## observation indicator
chn_all$obs <- 1:nrow(chn_all)

## interaction variables
chn_all$agereg2 <- interaction(chn_all$agegp_name,chn_all$reg2)
chn_all$agecause <- interaction(chn_all$agegp_name,chn_all$cause_name)
chn_all$reg2cause <- interaction(chn_all$reg2,chn_all$cause_name)
chn_all$agereg2cause <- interaction(chn_all$agegp_name,chn_all$reg2,chn_all$cause_name)
chn_all$rw_index <- as.numeric(interaction(chn_all$ageRW,chn_all$reg2,chn_all$causeRW))

## RW shared SDs
chn_all$rw_sd_index <- NA
chn_all$rw_sd_index[chn_all$cause_name=="injuries" & chn_all$agegp_name == "0-6d"] <- 1
chn_all$rw_sd_index[is.na(chn_all$rw_sd_index) & chn_all$reg2=="east.urban"] <- 2
chn_all$rw_sd_index[is.na(chn_all$rw_sd_index) & (chn_all$reg2=="west.rural" | chn_all$reg2=="mid.rural") & chn_all$agegp_name == "0-6d"] <- 3
chn_all$rw_sd_index[is.na(chn_all$rw_sd_index) & (chn_all$reg2=="west.rural" | chn_all$reg2=="mid.rural") & !(chn_all$agegp_name == "0-6d")] <- 4
chn_all$rw_sd_index[is.na(chn_all$rw_sd_index) & chn_all$reg2=="mid.urban" & chn_all$agegp_name == "0-6d"] <- 5
chn_all$rw_sd_index[is.na(chn_all$rw_sd_index)] <- chn_all$rw_index[is.na(chn_all$rw_sd_index)]
chn_all$rw_sd_index <- as.numeric(factor(chn_all$rw_sd_index))

##########################
## Fit model in INLA
##########################

# prior distributions (don't use these if you wish to fit the default priors in INLA)
rw_prec_prior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
iid_prec_prior <- list(prec = list(prior = "pc.prec", param = c(5, 0.01)))

chn_all <- chn_all[order(chn_all$rw_index,chn_all$agereg2cause,chn_all$year),]
mod_same_sd <- inla(deaths ~ factor(agegp)*factor(reg2) +
                        factor(agegp)*factor(cause) +
                        factor(reg2)*factor(cause) + 
                        f(year, replicate=rw_index, model = "rw2", constr = TRUE,
                          hyper = rw_prec_prior) +
                        f(factor(obs), model = "iid",
                          hyper = iid_prec_prior) +
                        offset(logpy),
                    family="poisson", 
                    data=chn_all,
                    quantiles = quantiles,
                    control.predictor=list(compute=TRUE),
                    control.compute=list(config = TRUE))

## save model results
setwd(savedir)
saveRDS(mod_same_sd,"results/china_results_pcpriors.RDS")

# load data (if only plotting)
mod_same_sd <- readRDS("results/china_results_pcpriors.RDS")

## load old data
mod_same_sd_old <- readRDS("results/china_results.RDS")

## get preds from old data
samples_old <- inla.posterior.sample(nsamples,mod_same_sd_old)
contents_old <- mod_same_sd_old$misc$configs$contents

id.pred_old <- which(contents_old$tag=="Predictor")
ind.pred_old <- contents_old$start[id.pred_old]-1 + (1:contents_old$length[id.pred_old])

samples.pred_old <- lapply(samples_old, function(x) {
    x$latent[ind.pred_old] - chn_all$logpy
})
s.pred_old <- t(matrix(unlist(samples.pred_old), byrow = T, nrow = length(samples.pred_old)))

##################
## Various posterior predictions
##################

## predictions for log mortaity rate
samples <- inla.posterior.sample(nsamples,mod_same_sd)

contents <- mod_same_sd$misc$configs$contents

id.pred <- which(contents$tag=="Predictor")
ind.pred <- contents$start[id.pred]-1 + (1:contents$length[id.pred])

samples.pred <- lapply(samples, function(x) {
    x$latent[ind.pred] - chn_all$logpy
    })
s.pred <- t(matrix(unlist(samples.pred), byrow = T, nrow = length(samples.pred)))
samples.pred <- NULL

## FE only predictions
id.fe <- which(contents$tag=="(Intercept)")
ind.fe <- contents$start[id.fe]-1 + (1:sum(contents$length[id.fe:length(contents$length)]))

X <- model.matrix(~factor(agegp)*factor(reg2) +
                      factor(agegp)*factor(cause) +
                      factor(reg2)*factor(cause), 
                  data = chn_all)
samples.fe <- lapply(samples, function(x) {
    as.vector(X %*% x$latent[ind.fe])
})
s.fe <- t(matrix(unlist(samples.fe), byrow = T, nrow = length(samples.fe)))
samples.fe <- NULL

## RW predictions
id.rw <- which(contents$tag=="year")
ind.rw <- contents$start[id.rw]-1 + (1:contents$length[id.rw])

samples.rw <- lapply(samples, function(x) {
    x$latent[ind.rw]
})
s.rw <- t(matrix(unlist(samples.rw), byrow = T, nrow = length(samples.rw)))
samples.rw <- NULL

## RE predictions
id.re <- which(contents$tag=="factor(obs)")
ind.re <- contents$start[id.re]-1 + (1:contents$length[id.re])

samples.re <- lapply(samples, function(x) {
    x$latent[ind.re]
})
s.re <- t(matrix(unlist(samples.re), byrow = T, nrow = length(samples.re)))
samples.re <- NULL

## FE + RW only predictions
s.fe.rw <- matrix(NA,nrow=nrow(chn_all),ncol=nsamples)
for (i in 1:nsamples) {
    if (i == 1) cat(paste("Starting sim 1 ...")); flush.console()
    if (i %% 10 == 0) cat(paste(i, "...")); flush.console()
    rwcount <- 0
    for (j in 1:length(unique(chn_all$rw_index))) {
        for (k in 1:nyear) {
            rwcount <- rwcount + 1
            s.fe.rw[chn_all$rw_index==j & chn_all$year==years[k],i] <- 
                s.fe[chn_all$rw_index==j & chn_all$year==years[k],i] + 
                s.rw[rwcount,i]
        }
    }
}

## predictions with added Poisson variation for predictions
pois_draws <- 100
s.pred.pois <- matrix(NA,nrow=nrow(chn_all),ncol=nsamples*pois_draws)
for (i in 1:nsamples) {
    if (i == 1) cat(paste("Starting sim 1 ...")); flush.console()
    if (i %% 10 == 0) cat(paste(i, "...")); flush.console()
    for (j in 1:pois_draws) {
        s.pred.pois[,((i-1)*(pois_draws))+j] <- rpois(nrow(chn_all),
                                                      rep(exp(s.pred[,i] + chn_all$logpy)))/exp(chn_all$logpy)
    }
}

## make results data frame
results <- chn_all

# fe only preds
results$alpha_m <- apply(s.fe,1,median)
results$alpha_l <- apply(s.fe,1,quantile,0.1)
results$alpha_u <- apply(s.fe,1,quantile,0.9)

## fe + rw preds
results$fe_rw_m <- apply(s.fe.rw,1,median)
results$fe_rw_l <- apply(s.fe.rw,1,quantile,0.1)
results$fe_rw_u <- apply(s.fe.rw,1,quantile,0.9)

# pred log mortality rates
results$pred_m <- apply(s.pred,1,median)
results$pred_l <- apply(s.pred,1,quantile,0.1)
results$pred_u <- apply(s.pred,1,quantile,0.9)

# pred deaths with Poisson noise
results$pred_mx_pois_m <- apply(s.pred.pois,1,median)
results$pred_mx_pois_l <- apply(s.pred.pois,1,quantile,0.1)
results$pred_mx_pois_u <- apply(s.pred.pois,1,quantile,0.9)

# pred log mortality rates from old non-pc prior model
results$pred_old_m <- apply(s.pred_old,1,median)
results$pred_old_l <- apply(s.pred_old,1,quantile,0.1)
results$pred_old_u <- apply(s.pred_old,1,quantile,0.9)

## comparison of old non pc prior and new models
results$pred_difference_old_new_m <- results$pred_m - results$pred_old_m

#-------------
# format data for saving
#-------------

## age in days variable
results$agegp_days <- as.numeric(as.character(factor(results$agegp,labels=c(3,18,90,255,540,1260))))

## change level order for age
results$agegp_name <- factor(results$agegp_name, levels = c("0-6d","7-27d","1-5m","6-11m","12-23m","24-59m"))

# rename regions in order to plot nicer
results$reg2 <- gsub("\\."," ",results$reg2)

# change age group labels
results <- results[order(results$reg2,results$agegp,results$year,results$cause),]
results$agegp_factor <- factor(results$agegp, labels = as.character(unique(results$agegp_name)))

## save results
saveRDS(results,"results/china_results_wide.RDS")

## load results (if only plotting)
results <- readRDS("results/china_results_wide.RDS")

# make into a tibble for making long dataset for plotting
res <- as_tibble(results)

# we are pivoting long to get a tidy dataset for easier plotting,
# and also getting rid of some unnecessary columns
res <- res %>% 
    pivot_longer(cols  = c("deaths", "exposure_old", "mx", "deaths_raw",
                           "exposure", "exposure_diff", "logmx", "logpy",
                           "alpha_m", "alpha_l", "alpha_u", 
                           "fe_rw_m", "fe_rw_l", "fe_rw_u", 
                           "pred_m", "pred_l", "pred_u", 
                           "pred_mx_pois_m", "pred_mx_pois_l", "pred_mx_pois_u", 
                           "pred_old_m", "pred_old_l", "pred_old_u", 
                           "pred_difference_old_new_m"),
                 names_to = "measure",
                 values_to = "value") %>%
    select("year", "agegp", "reg2", "cause", 
           "agegp_name", "year_name", "cause_name", "cause_ind", 
           "ageRW", "causeRW", 
           "agereg2", "agecause", "reg2cause", "agereg2cause", 
           "rw_index", "rw_sd_index",
           "agegp_days",     
           "measure", "value")

# calculate cause fractions for plotting over time
# pred
res2 <- res %>% 
    filter(measure == "fe_rw_m") %>%
    group_by(year, agegp, reg2) %>%
    mutate(value = exp(value)/sum(exp(value)))
res2$measure <- "pred_csmf_m"
res <- bind_rows(res,res2)

# SD of IID RE
sigma_e <- sqrt(1/mod_same_sd$summary.hyperpar["Precision for factor(obs)","0.5quant"])

# calculate residuals
res <- res %>% 
    pivot_wider(names_from = measure, values_from = value) %>%
    mutate(resid_std_mu = (deaths - exp(fe_rw_m + logpy))/sqrt(exp(fe_rw_m + logpy))) %>%
    mutate(resid_std_mu_sigma_e = (deaths - exp(fe_rw_m + logpy))/sqrt((exp(fe_rw_m + logpy)*exp(sigma_e^2)))) %>%
    pivot_longer(cols = c("deaths", "exposure_old", "mx", "deaths_raw",
                          "exposure", "exposure_diff", "logmx", "logpy",
                          "alpha_m", "alpha_l", "alpha_u", 
                          "fe_rw_m", "fe_rw_l", "fe_rw_u", 
                          "pred_m", "pred_l", "pred_u", 
                          "pred_mx_pois_m", "pred_mx_pois_l", "pred_mx_pois_u", 
                          "pred_old_m", "pred_old_l", "pred_old_u", 
                          "pred_difference_old_new_m",
                          "pred_csmf_m",
                          "resid_std_mu","resid_std_mu_sigma_e"), 
                 names_to = "measure", values_to = "value")

# save long dataset for some plots
saveRDS(res,"results/china_results_long.RDS")

results <- as.data.frame(results)

##########
## plot all log mortality rate results
##########

## plotting function
plotfun_mort <- function(results,data_var,
                         cause,cause_name,
                         region,
                         strat_var="agegp",by_var="year",
                         pub=FALSE) {
    
    ## for testing
    # cause <- "nCH16"
    # cause_name <- "other non-communicable"
    # region <- "east urban"
    # strat_var <- "agegp"
    # by_var <- "year"
    # data_var <- "logmx"
    # pub <- FALSE
    
    ## plotting size parameters
    if (pub) {
        cex.main <- 2.1
        cex.lab <- 2.2
        cex.axis <- 1.9
        cex.ylab <- cex.axis
        cex.legend <- 1.8
        cex <- 1.5
        cex.pt <- 1.7
    } else {
        cex.main <- 1.5
        cex.lab <- 1.5
        cex.axis <- 1.5
        cex.ylab <- cex.axis
        cex.legend <- 1.5
        cex <- 1.25
        cex.pt <- 1.5
    }
    
    ## set colors
    mycols <- brewer.pal(4,"Dark2")
    
    ## set plotting window
    if (pub) margins <- c(3.5,4.5,2,1) else margins <- c(3.5,4.5,4,1)
    par(mfrow=c(ceiling(length(unique(results[,strat_var]))/ceiling(sqrt(length(unique(results[,strat_var]))))),
                ceiling(sqrt(length(unique(results[,strat_var]))))),
        mar = margins)
    i <- 0
    
    if (strat_var=="year") strat_var_names <- c(1996:2015) else strat_var_names <- c("0-6d","7-27d","1-5m","6-11m","12-23m","24-59m")
    
    for (sv in unique(results[,strat_var])) {
        
        ## testing
        # sv <- 1
        
        if (strat_var == "agegp" & sv %in% c(4,5,6) & cause %in% c("nCH10","nCH11")) next
        
        i <- i+1
        
        ## index for what to plot
        index <- results$reg2==region & results[,strat_var]==sv & results$cause==cause
        
        ## indicator of 0 mortality points for plotting
        pch_plot <- ifelse(results[index,data_var] == -Inf,0,19)
        
        ## so we can plot 0 mortality points
        tmp_results <- results
        if (max(tmp_results[index,data_var]) == -Inf) {
            tmp_results[index & tmp_results[,data_var] == -Inf,data_var] <- min(tmp_results[index,"pred_l"]) - 0.25
        } else {
            tmp_results[index & tmp_results[,data_var] == -Inf,data_var] <- min(c(tmp_results[index & tmp_results[,data_var] != -Inf,data_var],
                                                                                  tmp_results[index,"pred_l"])) - 0.25
        }
        
        ## plot title
        if (pub) {
            plot_title <- strat_var_names[which(unique(results[,strat_var])==sv)]
        } else {
            if (i==1) {
                plot_title <- paste(paste(region,cause_name,sep="; "),"\n",
                                    strat_var_names[which(unique(results[,strat_var])==sv)],
                                    sep="")  
            } else {
                plot_title <- paste("\n",
                                    strat_var_names[which(unique(results[,strat_var])==sv)],
                                    sep="")   
            }
        }
        
        ## x axis to plot at
        if (by_var == "agegp") by_var_at <- results[index,"agegp_days"] else if (by_var == "year") by_var_at <- results[index,"year"]
        
        ## labels
        if (data_var == "logmx") ylabel <- "log(mortality rate)" else ylabel <- "mortality rate"
        
        if (i %in% c(1,4)) ylabel <- ylabel else ylabel <- ""
        ## plot mortality across byvar
        plot(tmp_results[index,"pred_m"]~by_var_at,type="l",lwd=2,lty=2,
             ylim=range(c(results[index,"pred_l"],
                          results[index,"pred_u"],
                          results[index,"alpha_l"],
                          results[index,"alpha_u"],
                          results[index,"fe_rw_l"],
                          results[index,"fe_rw_u"],
                          tmp_results[index,data_var])),
             main=plot_title,
             xlab="",ylab=ylabel,
             col=mycols[2],
             xaxt="n",
             cex=cex,
             cex.main=cex.main,
             cex.axis=cex.axis,
             cex.lab=cex.lab)
        
        if (by_var=="year") { 
            axislabels <- c(as.character(96:99),paste0("0",0:9),as.character(10:15))
        } else {
            axislabels <- c("0-6d","7-27d","1-5m","6-11m","12-23m","24-59m")
        }
        axis(1,at=by_var_at,labels=NA)
        text(x=by_var_at, y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
             labels=axislabels[1:length(by_var_at)], srt=45, adj=1, xpd=TRUE,
             cex=cex.ylab)
        lines(results[index,"alpha_m"]~by_var_at,col=mycols[3],lwd=2,lty=3)
        polygon(c(by_var_at,rev(by_var_at)),c(results[index,"alpha_l"],rev(results[index,"alpha_u"])),col=alpha(mycols[3],0.3),border=NA)
        lines(results[index,"fe_rw_m"]~by_var_at,col=mycols[1],lwd=2,lty=1)
        polygon(c(by_var_at,rev(by_var_at)),c(results[index,"fe_rw_l"],rev(results[index,"fe_rw_u"])),col=alpha(mycols[1],0.3),border=NA)
        polygon(c(by_var_at,rev(by_var_at)),c(results[index,"pred_l"],rev(results[index,"pred_u"])),col=alpha(mycols[2],0.3),border=NA)
        points(tmp_results[index,data_var]~by_var_at,pch=pch_plot,col=alpha(mycols[4],0.75),cex=cex.pt)
        if (i == 1) {
            legend("topright",c("final estimates","with IID REs","FEs only","observed data"),
                   fill=rep(NA,4),border=c(NA,NA,NA,NA),lwd=c(2,2,2,NA),
                   pch=c(NA,NA,NA,19),lty=c(1,2,3,NA),col=mycols[1:4],cex=cex.legend,
                   y.intersp=0.9)
            legend("topright",c("","","",""),lwd=c(NA,NA,NA,NA),
                   col=c(alpha(mycols[1],0.3),alpha(mycols[2],0.3),alpha(mycols[3],0.3),NA),border=c(NA,NA,NA,NA),
                   pch = rep(15,4),pt.cex = rep(3.85, 4),lty=c(NA,NA,NA,NA),cex=cex.legend,
                   y.intersp=0.9,x.intersp=8.38,xjust = 0,bty="n")
        }
    }
}

for (pp in c(TRUE,FALSE)) {
    ## testing
    # pp <- TRUE
    
    if (pp) {
        cat(paste("PUBLICATION VERSION:\n\n")); flush.console()
    } else {
        cat(paste("NON-PUBLICATION VERSION:\n\n")); flush.console()
    } 

    ## start pdf
    pdf(paste0("graphs/pub_",pp,"_logmx_ci_chn_inla_pcpriors_",Sys.Date(),".pdf"), 
        width=16,height=9)
    
    ## loop through regions and causes
    for (r in unique(results$reg2)) {
        for (c in unique(results$cause)) {
            
            ## testing
            # r <- unique(results$reg2)[1]
            # c <- unique(results$cause)[1]
            
            # cause name
            name_of_cause <- as.character(results$cause_name[results$cause==c][1])
            
            # make plot
            cat(paste(r,c,"\n",sep="; ")); flush.console()
            plotfun_mort(results,"logmx",cause=c,cause_name=name_of_cause,region=r,strat_var="agegp",by_var="year",pub=pp)
        }
    }
    
    dev.off()
}

#################
## Predicted death counts with additional Poisson variation
#################

pdf(paste0("graphs/mx_pi_chn_inla_pcpriors_",Sys.Date(),".pdf"), 
    width=16,height=8)

for (r in unique(results$reg2)) {
    for (c in unique(results$cause_name)) {
        cat(paste(r,c,"\n",sep="; ")); flush.console()
        p <- ggplot() + 
            geom_line(data = results[results$reg2 == r & results$cause_name == c,], 
                  aes(x = year_name, y = pred_mx_pois_m),
                  color = "dodgerblue") + 
            geom_ribbon(data = results[results$reg2 == r & results$cause_name == c,],
                        aes(x = year_name, 
                            ymin = pred_mx_pois_l, 
                            ymax = pred_mx_pois_u),
                        fill = "dodgerblue",
                        alpha = 0.2) +
            geom_point(data = results[results$reg2 == r & results$cause_name == c,],
                       aes(x = year_name, y = mx),
                       color = "tomato") +
            ylab("Mortality rate") + 
            xlab("Year") + 
            ggtitle(paste(r, c, sep = "; ")) +
            facet_wrap(~ agegp_name, scales = 'free') +
            theme_light()
        print(p)
    }
}

dev.off()

###############
## plot non-pc prior estimates
###############

## plotting function
plotfun_mort_old <- function(results,data_var,
                             cause,cause_name,
                             region,
                             strat_var="agegp",by_var="year",
                             pub=FALSE) {
    
    ## for testing
    # cause <- "nCH16"
    # cause_name <- "other non-communicable"
    # region <- "east urban"
    # strat_var <- "agegp"
    # by_var <- "year"
    # data_var <- "logmx"
    
    ## plotting size parameters
    if (pub) {
        cex.main <- 2.1
        cex.lab <- 2.2
        cex.axis <- 1.9
        cex.ylab <- cex.axis
        cex.legend <- 1.8
        cex <- 1.5
        cex.pt <- 1.7
    } else {
        cex.main <- 1.5
        cex.lab <- 1.5
        cex.axis <- 1.5
        cex.ylab <- cex.axis
        cex.legend <- 1.5
        cex <- 1.25
        cex.pt <- 1.5
    }
    
    ## set colors
    mycols <- brewer.pal(4,"Dark2")
    
    ## set plotting window
    if (pub) margins <- c(3.5,4.5,2,1) else margins <- c(5,4,4,1)+0.1
    par(mfrow=c(ceiling(length(unique(results[,strat_var]))/ceiling(sqrt(length(unique(results[,strat_var]))))),
                ceiling(sqrt(length(unique(results[,strat_var]))))),
        mar = margins)
    i <- 0
    
    if (strat_var=="year") strat_var_names <- c(1996:2015) else strat_var_names <- c("0-6d","7-27d","1-5m","6-11m","12-23m","24-59m")
    
    for (sv in unique(results[,strat_var])) {
        
        ## testing
        # sv <- 1
        
        if (strat_var == "agegp" & sv %in% c(4,5,6) & cause %in% c("nCH10","nCH11")) next
        
        i <- i+1
        
        ## index for what to plot
        index <- results$reg2==region & results[,strat_var]==sv & results$cause==cause
        
        ## indicator of 0 mortality points for plotting
        pch_plot <- ifelse(results[index,data_var] == -Inf,0,19)
        
        ## so we can plot 0 mortality points
        tmp_results <- results
        if (max(tmp_results[index,data_var]) == -Inf) {
            tmp_results[index & tmp_results[,data_var] == -Inf,data_var] <- min(tmp_results[index,"pred_l"]) - 0.25
        } else {
            tmp_results[index & tmp_results[,data_var] == -Inf,data_var] <- min(c(tmp_results[index & tmp_results[,data_var] != -Inf,data_var],
                                                                                  tmp_results[index,"pred_l"])) - 0.25
        }
        
        ## plot title
        if (pub) {
            plot_title <- strat_var_names[which(unique(results[,strat_var])==sv)]
        } else {
            if (i==1) {
                plot_title <- paste(paste(region,cause_name,sep="; "),"\n",
                                    strat_var_names[which(unique(results[,strat_var])==sv)],
                                    sep="")  
            } else {
                plot_title <- paste("\n",
                                    strat_var_names[which(unique(results[,strat_var])==sv)],
                                    sep="")   
            }
        }
        
        ## x axis to plot at
        if (by_var == "agegp") by_var_at <- results[index,"agegp_days"] else if (by_var == "year") by_var_at <- results[index,"year"]
        
        ## labels
        if (data_var == "logmx") ylabel <- "log(mortality rate)" else ylabel <- "mortality rate"
        
        if (i %in% c(1,4)) ylabel <- ylabel else ylabel <- ""
        ## plot mortality across byvar
        plot(tmp_results[index,"pred_old_m"]~by_var_at,type="l",lwd=2,
             ylim=range(c(results[index,"pred_old_l"],
                          results[index,"pred_old_u"],
                          tmp_results[index,data_var])),
             main=plot_title,
             xlab="",ylab=ylabel,
             col=mycols[1],
             xaxt="n",
             cex=cex,
             cex.main=cex.main,
             cex.axis=cex.axis,
             cex.lab=cex.lab)
        
        if (by_var=="year") { 
            axislabels <- c(as.character(96:99),paste0("0",0:9),as.character(10:15))
        } else {
            axislabels <- c("0-6d","7-27d","1-5m","6-11m","12-23m","24-59m")
        }
        axis(1,at=by_var_at,labels=NA)
        text(x=by_var_at, y=par()$usr[3]-0.04*(par()$usr[4]-par()$usr[3]),
             labels=axislabels[1:length(by_var_at)], srt=45, adj=1, xpd=TRUE,
             cex=cex.ylab)
        polygon(c(by_var_at,rev(by_var_at)),c(results[index,"pred_old_l"],rev(results[index,"pred_old_u"])),col=alpha(mycols[1],0.3),border=NA)
        points(tmp_results[index,data_var]~by_var_at,pch=pch_plot,col=alpha(mycols[4],0.75),cex=cex.pt)
        if (i == 1) {
            legend("topright",c("Final estimates","observed data"),
                   fill=rep(NA,2),border=c(NA,NA),lwd=c(2,NA),
                   pch=c(NA,19),lty=c(1,NA),col=c(mycols[1],mycols[4]),cex=cex.legend,
                   y.intersp=0.9)
            legend("topright",c("",""),lwd=c(NA,NA),
                   col=c(alpha(mycols[1],0.3),NA),border=c(NA,NA),
                   pch = rep(15,2),pt.cex = rep(4, 2),lty=c(NA,NA),cex=cex.legend,
                   y.intersp=0.9,x.intersp=8.42,xjust = 0,bty="n")
        }
    }
}

for (pp in c(TRUE,FALSE)) {
    ## testing
    # pp <- TRUE
    
    if (pp) {
        cat(paste("PUBLICATION VERSION:\n\n")); flush.console()
    } else {
        cat(paste("NON-PUBLICATION VERSION:\n\n")); flush.console()
    } 
    

    ## start pdf
    pdf(paste0("graphs/pub_",pp,"_logmx_ci_chn_inla_defaultpriors_",Sys.Date(),".pdf"), 
        width=16,height=9)
    
    ## loop through regions and causes
    for (r in unique(results$reg2)) {
        for (c in unique(results$cause)) {
            
            ## testing
            # r <- unique(results$reg2)[1]
            # c <- unique(results$cause)[1]
            
            # cause name
            name_of_cause <- as.character(results$cause_name[results$cause==c][1])
            
            # make plot
            cat(paste(r,c,"\n",sep="; ")); flush.console()
            plotfun_mort_old(results,"logmx",cause=c,cause_name=name_of_cause,region=r,
                             strat_var="agegp",by_var="year",pub=pp)
        }
    }
    
    dev.off()

}
