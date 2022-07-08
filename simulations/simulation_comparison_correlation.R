## Austin Schumacher
## 1/3/2019
## Simulation to motivate our method
## Goals
## - demonstrate problems with two-stage estimation framework
## - explore how a simple two-stage framework breaks in certain situations
## - compare this with a unified framework
## This code:
## - simulates data with fixed effects and correlated random effects
## - fits multistage and unified models
## - saves results

rm(list=ls())

## load libraries
library(mvtnorm); library(MCMCglmm); library(MASS);
library(foreign); library(gtools);library(boot);library(parallel);
library(scales);library(RColorBrewer);library(MGLM);
library(tidyverse); library(INLA);

## define directories

# directory for results
savedir <- "../../../Dropbox/SRS-child-mortality-output/"

# create folders to store results if necessary
if (!file.exists(paste0(savedir, "results"))) {
    dir.create(paste0(savedir, "results"))
}
if (!file.exists(paste0(savedir, "results/sims"))) {
    dir.create(paste0(savedir, "results/sims"))
}
if (!file.exists(paste0(savedir, "results/sims/tmp"))) {
    dir.create(paste0(savedir, "results/sims/tmp"))
}
if (!file.exists(paste0(savedir, "results/sims/tmp/correlation"))) {
    dir.create(paste0(savedir, "results/sims/tmp/correlation"))
}

## useful functions
expit <- function(x) exp(x)/(1+exp(x))

#######################
## SIMULATION SETTINGS
#######################

## code running
testing <- FALSE
casm_only <- FALSE

## which situation to simulate
intercept_only <- FALSE
linear_time <- FALSE
time_rw <- FALSE
interactions <- FALSE
multinomial <- FALSE

# load test data
chn_all <- readRDS("../data/mchss_test_data.RDS")

## format data
chn_all$logpy <- log(chn_all$exposure)
chn_all$ageRW <- ifelse(chn_all$agegp==1,1,2)
chn_all$causeRW <- as.character(chn_all$cause_name)
chn_all$causeRW[chn_all$cause_name=="diarrhea" | chn_all$cause_name=="other grp 1"] <- "diarrhea and other grp 1"
chn_all$causeRW[chn_all$cause_name=="congenital anomalies" | chn_all$cause_name=="other non-communicable"] <- "cong. anomalies and other n.c."

## data dimensions
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
nreg <- length(unique(chn_all$reg2))

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

## order data
chn_all <- chn_all[order(chn_all$rw_index,chn_all$agereg2cause,chn_all$year),]

## update data dimensions
ncause <- 2

## parameters, etc.
ndraws <- 1000 # number of draws for posterior estimates in all models

beta_intercepts <- c(-5,rep(0.5,nageg-1+nreg-1+ncause-1))
if (interactions) {
    tmp <- expand.grid(cause=1:ncause,ageg=1:nageg,reg=1:nreg,year=1:nyear)
    tmpX <- model.matrix(~ factor(cause) * factor(ageg) + factor(cause) * factor(reg),data=tmp)
    beta_intercepts <- c(-5,rep(0.5,ncol(tmpX)-1))
}
if (intercept_only) beta_intercepts <- c(-5,rep(0.5,ncause-1))
beta_time <- -0.05
gamma_sd <- 0.2
rw_sd <- 0.1
exposure <- median(chn_all$exposure)

nsim <- 100 # number of sims
rhos <- c(-0.5,0,0.5)
sigma_res <- c(0.01,0.1,1)

if (ncause > 5) stop("only 5 or fewer causes supported by INLA")

## start sim loop
for (sigma_re_ind in 1:length(sigma_res)) {
    
    if (testing) sigma_re_ind <- 1
    
    ## update status
    cat(paste0("\n Sigma: ",sigma_res[sigma_re_ind],"\n")); flush.console();
    
    # for (rho_ind in 1:1) {
    for (rho_ind in 1:length(rhos)) {
        
        if (testing) rho_ind <- 1
        
        ## update status
        cat(paste0("\n \t Correlation: ",rhos[rho_ind],"\n")); flush.console();
        
        sim_start <- 1
        n_alter <- 0
        
        for (sim in sim_start:(nsim+n_alter)) {
            
            if (testing) sim <- 88
            
            # adjust sim number for saving results due to skipping some sims taht error out
            sim_save <- sim - n_alter
            
            ## update status
            if (sim==1) cat(paste0("\n \t \t Starting sim 1 ...")); flush.console();
            if (sim %% 10 == 0) {
                cat(paste0(sim,"... ")); flush.console();
            }
            
            ## set seed
            set.seed(8008135+sim*2)
            
            ## parameters/indices
            sigma_re <- sigma_res[sigma_re_ind]
            rho <- rhos[rho_ind]
            
            ###################
            ## Simulate data
            ###################
            
            ## create data frame
            dat <- expand.grid(cause=1:ncause,ageg=1:nageg,reg=1:nreg,year=1:nyear)
            dat <- dat[order(dat$cause,dat$ageg,dat$reg,dat$year),]
            rownames(dat) <- NULL
            dat$intercept <- 1
            dat$year <- dat$year - mean(dat$year)
            years <- unique(dat$year)
            dat$logpy <- log(exposure)
            dat <- dat[order(dat$cause,dat$ageg,dat$reg,dat$year),]
            dat$index <- as.numeric(interaction(dat$year,dat$reg,dat$ageg))
            dat$obs <- 1:nrow(dat)
            nindex <- length(unique(dat$index))
            
            ## generate deaths
            fmla <- "~ factor(cause)"
            if (nageg > 1) fmla <- paste(fmla,"+ factor(ageg)")
            if (nreg > 1) fmla <- paste(fmla,"+ factor(reg)")
            if (linear_time) fmla <- paste(fmla,"+ year")
            if (interactions) {
                if (nageg > 1) {
                    fmla <- paste(fmla,"+ factor(ageg)*factor(cause)")
                }
                if (nreg > 1) {
                    fmla <- paste(fmla,"+ factor(reg)*factor(cause)")
                }
            }
            if (intercept_only) fmla <- "~ factor(cause)"
            fmla_mod <- paste("deaths",fmla)
            X <- model.matrix(as.formula(fmla), data=dat)
            
            # fixed effects
            beta <- beta_intercepts
            if (linear_time) beta <- c(beta_intercepts,beta_time)
            dat$alpha <- X%*%beta
            
            # correlated cause RE parameters
            if (ncause==2) {
                Omega <- matrix(c(1,rho,rho,1),2,2)
                true_rhos <- rho
                Sigma <- diag(rep(sigma_re,ncause),nrow = ncause,ncol = ncause) %*% 
                    Omega %*% 
                    diag(rep(sigma_re,ncause),nrow = ncause,ncol = ncause)
            } else {
                Omega <- matrix(0,nrow=ncause,ncol=ncause)
                diag(Omega) <- 1
                
                dcheck <- 1
                while (dcheck > 0) {
                    if (rho < 0) {
                        true_rhos <- runif(ncause*(ncause-1)/2,rho,0)
                    } else if (rho == 0) {
                        true_rhos <- rep(rho,ncause*(ncause-1)/2)
                    } else {
                        true_rhos <- runif(ncause*(ncause-1)/2,0,rho)
                    }
                    loopr <- 0
                    for (i in 1:ncause) {
                        for (j in 1:i) {
                            if (i==j) next
                            loopr <- loopr+1
                            Omega[i,j] <- Omega[j,i] <- true_rhos[loopr]
                        }
                    }
                    Sigma <- diag(rep(sigma_re,ncause),nrow = ncause,ncol = ncause) %*% 
                        Omega %*% 
                        diag(rep(sigma_re,ncause),nrow = ncause,ncol = ncause)
                    dcheck <- sum(eigen(Sigma)$values <= 0)
                }
            }
            
            true_rhos <- Omega[lower.tri(Omega)]
            
            ## generate data
            dat$gamma <- NA
            bmat <- rmvnorm(nindex,rep(0,ncause),Sigma)
            for (ii in 1:nindex) {
                for (cc in 1:ncause) {
                    dat$gamma[dat$index==ii & dat$cause == cc] <- bmat[ii,cc]
                }
            }
            
            # random walk for each cause
            Q <- matrix(NA,nrow=nyear,ncol=nyear)
            Q[1,] <- c(1,-2, 1,rep(0,nyear-3))
            Q[nyear,] <- c(rep(0,nyear-3),1,-2,1)
            Q[2,] <- c(-2,5,-4,1,rep(0,nyear-4))
            Q[nyear-1,] <- c(rep(0,nyear-4),1,-4,5,-2)
            for (i in 3:(nyear-2)) {
                Q[i,] <- c(rep(0,i-3),1,-4,6,-4,1,rep(0,nyear-2-i))
            }
            dat$rw <- NA
            for (rr in 1:ncause) {
                tau <- (rw_sd)^-2
                Qstar <- Q*tau
                eig <- eigen(Qstar)
                eval <- eig$values[c(-(nyear-1),-nyear)]
                evec <- eig$vectors[,c(-(nyear-1),-nyear)]
                randomwalk <- as.vector(evec%*%matrix(rnorm(length(eval),rep(0,length(eval)),1/sqrt(eval)),
                                                      nrow=length(eval),ncol=1))
                for (yy in 1:length(years)) {
                    dat$rw[dat$cause==rr & dat$year == years[yy]] <- randomwalk[yy]
                }
            }
            
            dat$logmx <- dat$alpha + dat$gamma
            if (time_rw) dat$logmx <- dat$logmx + dat$rw
            dat$deaths <- rpois(nrow(dat),exp(dat$logpy + dat$logmx))
            simdeathsallcause <- aggregate(dat$deaths,
                                           by=list(index=dat$index,
                                                   ageg=dat$ageg,
                                                   reg=dat$reg,
                                                   year=dat$year,
                                                   logpy=dat$logpy),
                                           sum)
            simdeathsallcause$deaths <- simdeathsallcause$x
            simdeathsallcause$x <- NULL
            simdeathsallcause <- simdeathsallcause[order(simdeathsallcause$index),]
            simdeathsallcause$obs <- 1:nrow(simdeathsallcause)
            
            # we get an error if all causes for an observation have 0 deaths
            zerodeaths <- unique(simdeathsallcause$index[simdeathsallcause$deaths == 0])
            simdeathsallcause <- simdeathsallcause[!(simdeathsallcause$index %in% zerodeaths),]
            dat <- dat[!(dat$index %in% zerodeaths),]
            nindex <- length(unique(dat$index))
            X_new <- model.matrix(as.formula(fmla), data=dat)
            
            ## store results
            results <- array(data=NA,dim=c(1,4))
            results_casm <- array(data=NA,dim=c(1,4))
            
            ## results to save from each exposure-simulation:
            ##  1. average relative bias for CSMRs
            ##  2. average absolute bias for CSMR
            ##  3. average coverage CSMRs
            ##  4. average width of CI for CSMR
            
            if (!casm_only) {
                
                #################
                ## Fit all-cause model
                #################
                
                ## all-cause model formula
                allcause_fmla <- gsub("~ factor\\(cause\\) \\+ ","~ ",fmla)
                allcause_fmla <- gsub("\\+ factor\\(ageg\\)\\*factor\\(cause\\)","",allcause_fmla)
                allcause_fmla <- gsub("\\+ factor\\(reg\\)\\*factor\\(cause\\)","",allcause_fmla)
                allcause_fmla <- gsub("factor\\(cause\\) \\+ ","",allcause_fmla)
                if (!interactions) allcause_fmla <- gsub("factor\\(cause\\)","1",allcause_fmla)
                
                X_allcause <- model.matrix(as.formula(allcause_fmla),data=simdeathsallcause)
                allcause_fmla_mod <- paste("deaths",allcause_fmla)
                allcause_fmla_mod <- paste(allcause_fmla_mod,"+ f(factor(obs), model='iid')")
                if (time_rw) allcause_fmla_mod <- paste(allcause_fmla_mod,"+ f(year,model='rw2',constr=TRUE)")
                allcause_fmla_mod <- paste0(allcause_fmla_mod," + offset(logpy)")
                
                # fit model
                allmod <- inla(as.formula(allcause_fmla_mod),data=simdeathsallcause,
                               family="poisson",
                               control.predictor=list(compute=TRUE),
                               control.compute=list(config=TRUE))
                
                ## sample from the posterior and get predicted inference
                samples <- inla.posterior.sample(ndraws,result=allmod,seed=nsim+sim*ndraws)
                
                contents <- allmod$misc$configs$contents
                
                id.pred <- which(contents$tag=="Predictor")
                ind.pred <- contents$start[id.pred]-1 + (1:contents$length[id.pred])
                
                samples.pred <- lapply(samples, function(x) {x$latent[ind.pred] - simdeathsallcause$logpy})
                
                s.pred <- matrix(unlist(samples.pred), byrow = T, nrow = length(samples.pred))
                colnames(s.pred) <- rownames(samples[[1]]$latent)[ind.pred]
                
                allcauseest <- t(s.pred)
                
                ##################
                ## Fit cause-specific model
                ##################
                
                if (multinomial) {
                    
                    csmf_fmla_mod <- gsub(" \\+ offset\\(logpy\\)","",allcause_fmla_mod)
                    csmf_fmla_mod <- gsub(" \\+ f\\(factor\\(obs\\), model='iid'\\) ","",csmf_fmla_mod)
                    csmf_fmla_mod <- paste0(csmf_fmla_mod," + f(ncause, model='iid', ",
                                            "constr = TRUE, hyper = ",
                                            "list(prec = list(initial = log(0.001),fixed = TRUE)))")
                    csmf_fmla_mod <- paste0(csmf_fmla_mod," + f(obs, model='iid",ncause,"d'",
                                            ", n = ncause * nindex")
                    
                    mod.csmf <- inla(as.formula(csmf_fmla_mod),
                                     family = "poisson",
                                     data = dat,
                                     control.predictor=list(compute=T),
                                     control.compute = list(config=TRUE))
                    
                    samples.csmf <- inla.posterior.sample(ndraws,result=mod.csmf,seed=nsim+sim*ndraws)
                    
                    contents.csmf <- mod.csmf$misc$configs$contents
                    
                    id.pred.csmf <- which(contents.csmf$tag=="Predictor")
                    ind.pred.csmf <- contents.csmf$start[id.pred.csmf]-1 + (1:contents.csmf$length[id.pred.csmf])
                    
                    samples.pred.csmf <- lapply(samples.csmf, 
                                                function(x) {
                                                    tmp <- data.frame(index=dat$index, 
                                                                      cause=dat$cause, 
                                                                      eta=exp(x$latent[ind.pred.csmf]))
                                                    (tmp %>% group_by(index) %>% mutate(prob = eta/sum(eta)))$prob
                                                })
                    
                    s.pred.csmf <- matrix(unlist(samples.pred.csmf), byrow = T, nrow = length(samples.pred.csmf))
                    colnames(s.pred.csmf) <- rownames(samples.csmf[[1]]$latent)[ind.pred.csmf]
                    
                    csmf_draws <- t(s.pred.csmf)
                } else {
                    tmpest <- c()
                    tmptotal <- matrix(0,ncol=ndraws,nrow=nindex)
                    for (i in 1:ncause) {
                        tmpcause <- dat[dat$cause==i,]
                        tmpmod <- inla(as.formula(allcause_fmla_mod),data=tmpcause,
                                       family="poisson",
                                       control.predictor=list(compute=TRUE),
                                       control.compute=list(config=TRUE))
                        
                        ## sample from the posterior and get predicted inference
                        samplestmp <- inla.posterior.sample(ndraws,result=tmpmod,seed=nsim+sim*ndraws)
                        
                        contentstmp <- tmpmod$misc$configs$contents
                        
                        id.pred.tmp <- which(contentstmp$tag=="Predictor")
                        ind.pred.tmp <- contentstmp$start[id.pred.tmp]-1 + (1:contentstmp$length[id.pred.tmp])
                        
                        samples.pred.tmp <- lapply(samplestmp, function(x) {exp(x$latent[ind.pred.tmp] - tmpcause$logpy)})
                        
                        s.pred.tmp <- matrix(unlist(samples.pred.tmp), byrow = T, nrow = length(samples.pred.tmp))
                        colnames(s.pred.tmp) <- rownames(samplestmp[[1]]$latent)[ind.pred.tmp]
                        
                        tmpest <- rbind(tmpest,t(s.pred.tmp))
                        tmptotal <- tmptotal + t(s.pred.tmp)
                    }
                    total <- tmptotal
                    for (i in 2:ncause) {
                        total <- rbind(total,tmptotal)
                    }
                    csmf_draws <- tmpest / total
                }
                
                ##################
                ## Combine all-cause and cause-specific draws
                ##################
                
                ## CSMR draws
                logcsmr_draws <- matrix(NA,nrow=nrow(dat),ncol=ndraws)
                allcause_draws <- exp(allcauseest)
                
                for (cc in 2:ncause) {
                    allcause_draws <- rbind(allcause_draws,exp(allcauseest))
                }
                logcsmr_draws <- log(allcause_draws*csmf_draws)
                
                ## inference
                # CIs
                csmr_est <- apply(logcsmr_draws,1,median)
                cis_csmr <- apply(logcsmr_draws,1,quantile,c(0.025,0.975))
                
                # CSMR results
                results[,1] <- mean((csmr_est- dat$logmx)/dat$logmx)
                results[,2] <- mean((csmr_est- dat$logmx))
                results[,3] <- mean(dat$logmx > cis_csmr[1,] & dat$logmx < cis_csmr[2,])
                results[,4] <- mean(cis_csmr[2,] - cis_csmr[1,])
            }
            
            ##################
            ## Unified model
            ##################
            
            fmla.csmr.mod <- fmla_mod
            fmla.csmr.mod <- paste0(fmla.csmr.mod," + f(obs, model='iid",ncause,"d', n=ncause*nindex)")
            if (time_rw) fmla.csmr.mod <- paste(fmla.csmr.mod,"+ f(year,model='rw2',constr=TRUE)")
            fmla.csmr.mod <- paste(fmla.csmr.mod,"+ offset(logpy)")
            
            fit1 <- inla(as.formula(fmla.csmr.mod),
                         data=dat,
                         family="poisson",
                         control.predictor=list(compute=TRUE),
                         control.compute=list(config=TRUE))
            pred_logmx <- log(fit1$summary.fitted.values$`0.5quant`/exp(dat$logpy))
            cis_logmx <- cbind(log(fit1$summary.fitted.values$`0.025quant`/exp(dat$logpy)),
                               log(fit1$summary.fitted.values$`0.975quant`/exp(dat$logpy)))
            
            # CSMR
            results_casm[,1] <- mean((pred_logmx - dat$logmx)/dat$logmx)
            results_casm[,2] <- mean((pred_logmx - dat$logmx))
            results_casm[,3] <- mean((dat$logmx > cis_logmx[,1]) & (dat$logmx < cis_logmx[,2]))
            results_casm[,4] <- mean(cis_logmx[,2] - cis_logmx[,1])
            
            ################
            ## Save results
            ################
            
            if (!testing) {
                ## save results
                if (!casm_only) {
                    saveRDS(results,
                            file=paste0(savedir, 
                                        "results/sims/tmp/correlation/motive_corr_res",
                                        "_intonly_",intercept_only,
                                        "_lineartime_",linear_time,
                                        "_timerw_",time_rw,
                                        "_interactions_",interactions,
                                        "_multinom_",multinomial,
                                        "_ncause_",ncause,
                                        "_sim_",sim_save,
                                        "_sigma_",sigma_re,
                                        "_rho_",rho,
                                        ".rds"))
                }
                saveRDS(results_casm,
                        file=paste0(savedir,
                                    "results/sims/tmp/correlation/motive_corr_res_casm",
                                    "_intonly_",intercept_only,
                                    "_lineartime_",linear_time,
                                    "_timerw_",time_rw,
                                    "_interactions_",interactions,
                                    "_multinom_",multinomial,
                                    "_ncause_",ncause,
                                    "_sim_",sim_save,
                                    "_sigma_",sigma_re,
                                    "_rho_",rho,
                                    ".rds"))
            }
        }
    }
}
