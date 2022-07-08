# Austin Schumacher
# 11/20/2019
# This code:
# - simulates data with correlated random effects
# - fits unified models (in Stan) to the data
# - plots simulation parameters vs. the estimates

rm(list=ls())

## load libraries
library(Rcpp); library(StanHeaders); library(BH); library(rstan);library(bayesplot);
library(mvtnorm); library(MASS);library(gtools); library(parallel);
library(scales); library(RColorBrewer);library(data.table);
library(ggplot2);

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
if (!file.exists(paste0(savedir, "results/sims/explore_correlation"))) {
    dir.create(paste0(savedir, "results/sims/explore_correlation"))
}
if (!file.exists(paste0(savedir, "graphs"))) {
    dir.create(paste0(savedir, "graphs"))
}
if (!file.exists(paste0(savedir, "graphs/sims"))) {
    dir.create(paste0(savedir, "graphs/sims"))
}
if (!file.exists(paste0(savedir, "graphs/sims/explore_correlation"))) {
    dir.create(paste0(savedir, "graphs/sims/explore_correlation"))
}

#######################
## SETTINGS
#######################

## code running options
testing <- FALSE

## data generation options
cause_re <- FALSE
few_corrs <- FALSE
time_rw <- TRUE
no_corr <- FALSE
FE_interactions <- TRUE

## modeling options
cause_re_model <- FALSE
time_rw_model <- TRUE
no_corr_model <- FALSE # IMPLIES CAUSE RE AND TIME RW BUT NO CORRELATIONS
FE_interactions_model <- TRUE

## set seed
set.seed(881352)

## set number of cores for parallel computing
ncores <- detectCores()

## parameters
nyears <- c(20)
nreg <- 6
nageg <- 6
ncause <- 8
sigma_rws <- c(0.1)
sigma_res <- c(0.5)
corrs <- c(0.75)
quantiles <- c(0.1,0.5,0.9)
rw_vars <- c("ageg","cause")
re_vars <- c("ageg","reg","cause")
exposure <- 10000
lkj_hyperpar <- 1

## STAN options
rstan_options(auto_write = TRUE)
niter <- 4000
nchains <- 2
nthin <- 1
prop_warmup <- 1/2
if (testing) niter <- 200
if (testing) nchains <- 2
if (testing) nthin <- 1
if (testing) prop_warmup <- 5/10

for (yy in 1:length(nyears)) {
    for (ss in 1:length(sigma_rws)) {
        for (ssre in 1:length(sigma_res)) {
            for (cc in 1:length(corrs)) {
                
                if (testing) {
                    yy <- ss <- 1
                    ssre <- 1
                    cc <- 1
                }
                
                nyear <- nyears[yy]
                sigma_rw <- sigma_rws[ss]
                sigma_re <- sigma_res[ssre]
                rho <- corrs[cc]
                # sigma_epsilon <- sigma_epsilons[sse]
                tau_rw <- 1/(sigma_rw^2)
                # tau_epsilon <- 1/(sigma_epsilon^2)
                
                ## create data frame
                dat <- expand.grid(reg=1:nreg,cause=1:ncause,ageg=1:nageg,year=1:nyear)
                dat <- dat[order(dat$cause,dat$ageg,dat$reg,dat$year),]
                rownames(dat) <- NULL
                dat$intercept <- 1
                dat$re_index <- as.numeric(interaction(dat[,re_vars]))
                n_re_index <- length(unique(dat$re_index))
                dat$rw_index <- as.numeric(interaction(dat[,rw_vars]))
                n_rw_index <- length(unique(dat$rw_index))
                dat$logpy <- log(exposure)
                
                ## generate RW
                if (time_rw) {
                    Q <- matrix(NA,nrow=nyear,ncol=nyear)
                    Q[1,] <- c(1,-2, 1,rep(0,nyear-3))
                    Q[nyear,] <- c(rep(0,nyear-3),1,-2,1)
                    Q[2,] <- c(-2,5,-4,1,rep(0,nyear-4))
                    Q[nyear-1,] <- c(rep(0,nyear-4),1,-4,5,-2)
                    for (i in 3:(nyear-2)) {
                        Q[i,] <- c(rep(0,i-3),1,-4,6,-4,1,rep(0,nyear-2-i))
                    }
                    Q <- Q*tau_rw
                    eig <- eigen(Q)
                    eval <- eig$values[c(-(nyear-1),-nyear)]
                    evec <- eig$vectors[,c(-(nyear-1),-nyear)]
                    dat$x_t <- NA
                    for (rr in 1:n_rw_index) {
                        randomwalk <- as.vector(evec%*%matrix(rnorm(length(eval),rep(0,length(eval)),1/sqrt(eval)),nrow=length(eval),ncol=1))
                        randomwalk <- randomwalk - mean(randomwalk)
                        for (yr in 1:nyear) {
                            dat$x_t[dat$rw_index==rr & dat$year==yr] <- randomwalk[yr]
                        }
                    }
                } else {
                    dat$x_t <- 0
                }
                
                ## correlated cause RE parameters
                if (cause_re) {
                    if (no_corr) {
                        Omega <- matrix(0,nrow=ncause,ncol=ncause)
                        diag(Omega) <- 1
                    } else {
                        if (ncause==2) {
                            Omega <- matrix(c(1,rho,rho,1),2,2)
                            true_rhos <- rho
                        } else {
                            Omega <- matrix(0,nrow=ncause,ncol=ncause)
                            diag(Omega) <- 1
                            
                            dcheck <- 1
                            while (dcheck >0) {
                                if (few_corrs) {
                                    if (ncause == 3 | ncause  == 4) {
                                        Omega[1,2] <- Omega[2,1] <- rho
                                    } else if (ncause == 5 | ncause == 6) {
                                        Omega[1,2] <- Omega[2,4] <- Omega[2,1] <- Omega[4,2] <- rho
                                        Omega[3,5] <- Omega[5,3] <- -1*rho
                                    } else {
                                        Omega[1,2] <- Omega[2,4] <- Omega[2,1] <- Omega[4,2] <- rho
                                        Omega[3,7] <- Omega[7,3] <- -1*rho
                                    }
                                    true_rhos <- as.vector(Omega[lower.tri(Omega)])
                                } else {
                                    if (rho < 0) {
                                        true_rhos <- runif(ncause*(ncause-1)/2,rho,0)
                                    } else if (rho == 0) {
                                        true_rhos <- rep(rho,ncause*(ncause-1)/2)
                                    } else {
                                        true_rhos <- runif(ncause*(ncause-1)/2,-1*rho,rho)
                                    }
                                    loopr <- 0
                                    for (i in 1:ncause) {
                                        for (j in 1:i) {
                                            if (i==j) next
                                            loopr <- loopr+1
                                            Omega[i,j] <- Omega[j,i] <- true_rhos[loopr]
                                        }
                                    }
                                }
                                Sigma <- diag(rep(sigma_re,ncause)) %*% Omega %*% diag(rep(sigma_re,ncause))
                                dcheck <- sum(eigen(Sigma)$values <= 0)
                            }
                        }
                    }
                    
                    true_rhos <- Omega[lower.tri(Omega)]
                    
                    if (no_corr) {
                        dat$b <- NA
                        for (ii in 1:n_re_index) dat$b[dat$re_index==ii] <- rnorm(1,0,sigma_re)
                    } else {
                        bmat <- rmvnorm(n_re_index,rep(0,ncause),Sigma)
                        for (ii in 1:n_re_index) {
                            for (cc in 1:ncause) {
                                dat$b[dat$re_index==ii & dat$cause == cc] <- bmat[ii,cc]
                            }
                        }
                    }
                } else {
                    dat$b <- 0
                }
                
                if (FE_interactions) {
                    if (nreg == 1) {
                        FE_formula <- as.formula(~ factor(ageg) * factor(cause))
                    } else {
                        FE_formula <- as.formula(~ factor(ageg)*factor(reg) + factor(ageg)*factor(cause) + factor(reg)*factor(cause))
                    }
                } else {
                    if (nreg == 1) {
                        FE_formula <- as.formula(~ factor(ageg) + factor(cause))
                    } else {
                        FE_formula <- as.formula(~ factor(ageg) + factor(reg) + factor(cause))
                    }
                }
                
                X_truth <- model.matrix(FE_formula,data=dat)
                beta <- c(-6,rep(0.5,ncol(X_truth)-1))
                
                ## generate data
                dat$alpha <- as.vector(X_truth%*%beta)
                dat$y <- rpois(nrow(dat),exp(dat$logpy + dat$alpha + dat$b + dat$x_t))
                dat$mx <- dat$y/exp(dat$logpy)
                dat$logmx <- log(dat$mx)
                dat$lambda <- dat$alpha + dat$b + dat$x_t
                
                ###########
                ## Stan
                ###########
                
                if (FE_interactions_model) {
                    if (nreg == 1) {
                        FE_formula_model <- as.formula(~ factor(ageg) * factor(cause))
                    } else {
                        FE_formula_model <- as.formula(~ factor(ageg)*factor(reg) + factor(ageg)*factor(cause) + factor(reg)*factor(cause))
                    }
                } else {
                    if (nreg == 1) {
                        FE_formula_model <- as.formula(~ factor(ageg) + factor(cause))
                    } else {
                        FE_formula_model <- as.formula(~ factor(ageg) + factor(reg) + factor(cause))
                    }
                }
                
                X <- model.matrix(FE_formula_model,data=dat)
                
                if (cause_re_model & time_rw_model) {
                    datlist <- list(N = nrow(dat),
                                    p = ncol(X),
                                    X = X,
                                    ncause = ncause,
                                    s = n_re_index,
                                    n_rw_index = n_rw_index,
                                    nyear=nyear,
                                    s_index = dat$re_index,
                                    c = dat$cause,
                                    k = dat$year,
                                    rw_index = dat$rw_index,
                                    logpy = dat$logpy,
                                    y = dat$y,
                                    lkj_hyperpar = lkj_hyperpar)
                    
                    ## initial values
                    estimates <- function(y,logpy,s,s_index,ncause,c,X,n_rw_index,nyear,
                                          perturb=chain_id>1) {
                        
                        ## testing
                        # s_index = dat$re_index
                        # s <- n_re_index
                        # c = dat$cause
                        # logpy = dat$logpy
                        # y = dat$y
                        # perturb = FALSE
                        
                        if(perturb) {
                            y <- y + as.integer(round(rnorm(length(y), 0, 20)))
                            y[y < 0] <- 0
                        }
                        
                        tmp <- y+1
                        # tmp[tmp==0] <- 1
                        beta <- solve(t(X) %*% X, t(X) %*% log(tmp/exp(logpy)), tol = 1e-12)
                        if (sum(dim(beta))==2) beta <- as.vector(beta)
                        alpha <- X %*% beta
                        
                        alpha_mat <- matrix(NA,nrow=s,ncol=ncause)
                        b_star <- matrix(NA,nrow=s,ncol=ncause)
                        for (i in 1:length(tmp)) {
                            alpha_mat[s_index[i],c[i]] <- alpha[i]
                            b_star[s_index[i],c[i]] <- log(tmp[i]) - logpy[i] - alpha[i]
                        }
                        
                        Lcorr <- t(chol(cor(b_star)+diag(0.0001,ncause)))
                        sigma <- sqrt(diag(var(b_star)))
                        
                        gamma_mat <- matrix(0,nrow=n_rw_index,ncol=nyear)
                        sigma_gamma <- 0.5
                        
                        return(list(Lcorr=Lcorr, sigma=sigma, b_star=b_star, beta=beta,
                                    gamma_mat=gamma_mat,sigma_gamma=sigma_gamma))
                    }
                    
                    inits <- function(chain_id){
                        values <- estimates(datlist$y, datlist$logpy, datlist$s,datlist$s_index, 
                                            datlist$ncause, datlist$c, datlist$X,
                                            datlist$n_rw_index,datlist$nyear,
                                            perturb = TRUE)
                        values[["beta"]] <- as.vector(values[["beta"]])
                        return(values)
                    }
                    
                    mod_file <- "stan_models/model_causeRE_timeRW2constr.stan"
                    
                } else if (cause_re_model & !time_rw_model) {
                    
                    datlist <- list(N = nrow(dat),
                                    p = ncol(X),
                                    X = X,
                                    ncause = ncause,
                                    s = n_re_index,
                                    s_index = dat$re_index,
                                    c = dat$cause,
                                    logpy = dat$logpy,
                                    y = dat$y,
                                    lkj_hyperpar = lkj_hyperpar)
                    
                    ## initial values
                    estimates <- function(y,logpy,s,s_index,ncause,c,X,
                                          perturb=chain_id>1) {
                        
                        ## testing
                        # s_index = dat$re_index
                        # s <- n_re_index
                        # c = dat$cause
                        # logpy = dat$logpy
                        # y = dat$y
                        # perturb = FALSE
                        
                        if(perturb) {
                            y <- y + as.integer(round(rnorm(length(y), 0, 20)))
                            y[y < 0] <- 0
                        }
                        
                        tmp <- y+1
                        # tmp[tmp==0] <- 1
                        beta <- solve(t(X) %*% X, t(X) %*% log(tmp/exp(logpy)), tol = 1e-12)
                        if (sum(dim(beta))==2) beta <- as.vector(beta)
                        alpha <- X %*% beta
                        
                        alpha_mat <- matrix(NA,nrow=s,ncol=ncause)
                        b_star <- matrix(NA,nrow=s,ncol=ncause)
                        for (i in 1:length(tmp)) {
                            alpha_mat[s_index[i],c[i]] <- alpha[i]
                            b_star[s_index[i],c[i]] <- log(tmp[i]) - logpy[i] - alpha[i]
                        }
                        
                        Lcorr <- t(chol(cor(b_star)+diag(0.0001,ncause)))
                        sigma <- sqrt(diag(var(b_star)))
                        
                        return(list(Lcorr=Lcorr, sigma=sigma, b_star=b_star, beta=beta))
                    }
                    
                    inits <- function(chain_id){
                        values <- estimates(datlist$y, datlist$logpy, datlist$s,datlist$s_index, 
                                            datlist$ncause, datlist$c, datlist$X,
                                            perturb = TRUE)
                        values[["beta"]] <- as.vector(values[["beta"]])
                        return(values)
                    }
                    
                    mod_file <- "stan_models/model_causeRE_only_all_age.stan"
                    
                } else if (!cause_re_model & time_rw_model) {
                    
                    datlist <- list(N = nrow(dat),
                                    p = ncol(X),
                                    X = X,
                                    n_rw_index = n_rw_index,
                                    nyear=nyear,
                                    k = dat$year,
                                    rw_index = dat$rw_index,
                                    logpy = dat$logpy,
                                    y = dat$y)
                    
                    ## initial values
                    estimates <- function(y,logpy,X,n_rw_index,nyear,
                                          perturb=chain_id>1) {
                        
                        ## testing
                        # s_index = dat$re_index
                        # s <- n_re_index
                        # c = dat$cause
                        # logpy = dat$logpy
                        # y = dat$y
                        # perturb = FALSE
                        
                        if(perturb) {
                            y <- y + as.integer(round(rnorm(length(y), 0, 20)))
                            y[y < 0] <- 0
                        }
                        
                        tmp <- y+1
                        # tmp[tmp==0] <- 1
                        beta <- solve(t(X) %*% X, t(X) %*% log(tmp/exp(logpy)), tol = 1e-12)
                        if (sum(dim(beta))==2) beta <- as.vector(beta)
                        
                        gamma_mat <- matrix(0,nrow=n_rw_index,ncol=nyear)
                        sigma_gamma <- 0.5
                        
                        return(list(beta=beta,gamma_mat=gamma_mat,sigma_gamma=sigma_gamma))
                    }
                    
                    inits <- function(chain_id){
                        values <- estimates(datlist$y, datlist$logpy, datlist$X,
                                            datlist$n_rw_index,datlist$nyear,
                                            perturb = TRUE)
                        values[["beta"]] <- as.vector(values[["beta"]])
                        return(values)
                    }
                    
                    mod_file <- "stan_models/model_timeRW2constr.stan"
                    
                } else if (no_corr_model) {
                    datlist <- list(N = nrow(dat),
                                    p = ncol(X),
                                    X = X,
                                    s = n_re_index,
                                    n_rw_index = n_rw_index,
                                    nyear=nyear,
                                    s_index = dat$re_index,
                                    k = dat$year,
                                    rw_index = dat$rw_index,
                                    logpy = dat$logpy,
                                    y = dat$y)
                    
                    ## initial values
                    estimates <- function(y,logpy,s,s_index,X,n_rw_index,nyear,
                                          perturb=chain_id>1) {
                        
                        ## testing
                        s_index = dat$re_index
                        s <- n_re_index
                        logpy = dat$logpy
                        y = dat$y
                        perturb = TRUE
                        
                        if(perturb) {
                            y <- y + as.integer(round(rnorm(length(y), 0, 20)))
                            y[y < 0] <- 0
                        }
                        
                        tmp <- y+1
                        # tmp[tmp==0] <- 1
                        beta <- solve(t(X) %*% X, t(X) %*% log(tmp/exp(logpy)), tol = 1e-12)
                        if (sum(dim(beta))==2) beta <- as.vector(beta)
                        alpha <- X %*% beta
                        
                        alpha_s <- rep(NA,s)
                        b_star <- rep(NA,s)
                        for (i in 1:length(tmp)) {
                            alpha_s[s_index[i]] <- alpha[i]
                            b_star[s_index[i]] <- log(tmp[i]) - logpy[i] - alpha[i]
                        }
                        
                        sigma <- sd(b_star)
                        
                        gamma_mat <- matrix(0,nrow=n_rw_index,ncol=nyear)
                        sigma_gamma <- 0.5
                        
                        return(list(sigma=sigma, b_star=b_star, beta=beta,
                                    gamma_mat=gamma_mat,sigma_gamma=sigma_gamma))
                    }
                    
                    inits <- function(chain_id){
                        values <- estimates(datlist$y, datlist$logpy, datlist$s,datlist$s_index, 
                                            datlist$X,
                                            datlist$n_rw_index,datlist$nyear,
                                            perturb = TRUE)
                        values[["beta"]] <- as.vector(values[["beta"]])
                        return(values)
                    }
                    
                    mod_file <- "stan_models/model_nocorrRE_timeRW2constr.stan"
                } else {
                    stop("no model coded up for this yet")
                }
                
                # set number of cores to the max possible or the number of chains, whichever is smaller
                options(mc.cores = parallel::detectCores())
                
                # fit the model
                setwd(wd)
                
                start.time <- proc.time() 
                mod_stan <- stan(file = mod_file,
                                 data = datlist,
                                 iter = niter, chains = nchains, thin=nthin,
                                 warmup = niter*prop_warmup,
                                 init = inits,
                                 control = list(adapt_delta = 0.95, max_treedepth = 15))
                stop.time <- proc.time()
                elapsed.time <- stop.time[3] - start.time[3]

                parlist <- c("beta","log_lik")
                if (cause_re_model) {
                    parlist <- c(parlist,"sigma","b_star","Omega")
                } 
                if (time_rw_model) {
                    parlist <- c(parlist,"sigma_gamma","gamma_mat")
                } 
                if (no_corr_model) {
                    parlist <- c(parlist,"sigma","b_star","sigma_gamma","gamma_mat")
                }
                res_stan <- extract(mod_stan,pars=parlist,permuted=FALSE,inc_warmup=FALSE)

                ## we are saving the model information in a csv to reference the models we save
                if (file.exists(paste0(savedir, "results/sims/explore_correlation/allmods.csv"))) {
                    # read in current models
                    allmods <- read.csv(paste0(savedir, "results/sims/explore_correlation/allmods.csv"),header=TRUE,stringsAsFactors = FALSE)
                    
                    # get model number
                    model_number <- max(allmods$model_number) + 1
                } else {
                    model_number <- 1
                }
                addrow <- data.frame(model_number=model_number,
                                     date=as.character(Sys.time()),
                                     runtime = elapsed.time,
                                     mod_file = mod_file,
                                     testing = testing,
                                     cause_re = cause_re,
                                     few_corrs = few_corrs,
                                     time_rw = time_rw,
                                     no_corr = no_corr,
                                     FE_interactions = FE_interactions,
                                     cause_re_model = cause_re_model,
                                     time_rw_model = time_rw_model,
                                     no_corr_model = no_corr_model,
                                     FE_interactions_model = FE_interactions_model,
                                     FE_formula = as.character(FE_formula)[2],
                                     FE_formula_model = as.character(FE_formula_model)[2],
                                     sigma_rw = sigma_rw,
                                     sigma_re = sigma_re,
                                     rho = rho,
                                     lkj_hyperpar = lkj_hyperpar,
                                     nyear=nyear,
                                     ncause = ncause,
                                     nreg = nreg,
                                     nageg = nageg,
                                     exposure = exposure,
                                     re_vars = paste(re_vars,collapse="_"),
                                     rw_vars = paste(rw_vars,collapse="_"),
                                     nchains = nchains,
                                     niter = niter,
                                     prop_warmup = prop_warmup,
                                     nthin = nthin)
                
                if (file.exists(paste0(savedir, "results/sims/explore_correlation/allmods.csv"))) {
                    allmods <- rbind(allmods,addrow)
                } else {
                    allmods <- addrow
                }
                
                # save updated models document
                write.csv(allmods,file=paste0(savedir, "results/sims/explore_correlation/allmods.csv"),row.names = FALSE)
                
                # save model chains
                saveRDS(res_stan,file=paste0(savedir, "results/sims/explore_correlation/model_",model_number,"_res.RDS"))
                
                # save parameters
                if (cause_re & !time_rw) {
                    save(dat,beta,sigma_re,rho,Omega,re_vars, 
                         file = paste0(savedir, "results/sims/explore_correlation/model_",model_number,"_pars_data.Rdata"))
                } else if (time_rw & !cause_re) {
                    save(dat,beta,sigma_rw,rw_vars, 
                         file = paste0(savedir, "results/sims/explore_correlation/model_",model_number,"_pars_data.Rdata"))
                } else if ((time_rw & cause_re)) {
                    save(dat,beta,sigma_re,rho,Omega,re_vars,sigma_rw,rw_vars,
                         file = paste0(savedir, "results/sims/explore_correlation/model_",model_number,"_pars_data.Rdata"))
                } else if (no_corr) {
                    save(dat,beta,sigma_re,re_vars,sigma_rw,rw_vars,
                         file = paste0(savedir, "results/sims/explore_correlation/model_",model_number,"_pars_data.Rdata"))
                }
                
                # save initial values
                mod_inits <- get_inits(mod_stan)
                saveRDS(mod_inits,file=paste0(savedir, "results/sims/explore_correlation/model_",model_number,"_inits.RDS"))
                
                # save summaries
                mod_summary <- summary(mod_stan,pars=parlist[!(parlist %in% c("gamma_mat","b_star"))],probs=c(0.1,0.5,0.9))$summary
                saveRDS(mod_summary,file=paste0(savedir, "results/sims/explore_correlation/model_",model_number,"_summaries.RDS"))

                # save diagnostics
                mod_diag <- data.frame(n_divergent_iters=sum(get_divergent_iterations(mod_stan)),
                                       n_max_tree_depth_iters=sum(get_max_treedepth_iterations(mod_stan)),
                                       n_low_BFMI_chains=sum(which(get_bfmi(mod_stan)<0.2)))
                saveRDS(mod_diag,file=paste0(savedir, "results/sims/explore_correlation/model_",model_number,"_diag.RDS"))
                
                # save graphs of the model diagnostics
                pdf(paste0(savedir, "graphs/sims/explore_correlation/model_",model_number,"_traceplots.pdf"),width=10,height=10)
                print(stan_trace(mod_stan,pars=c("beta")))
                if (cause_re_model) {
                    print(stan_trace(mod_stan,pars=c("sigma")))
                    print(stan_trace(mod_stan,pars=c("Omega")))
                }
                if (time_rw_model) {
                    print(stan_trace(mod_stan,pars=c("sigma_gamma")))
                }
                if (no_corr_model) {
                    print(stan_trace(mod_stan,pars=c("sigma")))
                    print(stan_trace(mod_stan,pars=c("sigma_gamma")))
                }
                dev.off()
                
                # more plots of model diagnostics
                pdf(paste0(savedir, "graphs/sims/explore_correlation/model_",model_number,"_traceplotsREsample.pdf"),width=4,height=4)
                
                if (time_rw_model | no_corr_model) {
                    for (gg in sample(1:n_rw_index,ceiling(n_rw_index/10),replace=FALSE)) {
                        for (yy in sample(1:nyear,ceiling(nyear/10),replace=FALSE)) {
                            print(stan_trace(mod_stan,pars=paste0("gamma_mat[",gg,",",yy,"]")))
                        }
                    }
                }
                if (cause_re_model) {
                    for (ee in sample(1:n_re_index,ceiling(n_re_index/10),replace=FALSE)) {
                        for (ccc in sample(1:ncause,ceiling(ncause/3),replace=FALSE)) {
                            print(stan_trace(mod_stan,pars=paste0("b_star[",ee,",",ccc,"]")))
                        }
                    }
                }
                if (no_corr_model) {
                    for (ee in sample(1:n_re_index,ceiling(n_re_index/10),replace=FALSE)) {
                        print(stan_trace(mod_stan,pars=paste0("b_star[",ee,"]")))
                    }
                }
                dev.off()
                
                ## stan diagnostic plots
                pdf(paste0(savedir, "graphs/sims/explore_correlation/model_",model_number,"_diagplots.pdf"),width=10,height=10)
                stan_rhat(mod_stan)
                stan_diag(mod_stan)
                stan_ess(mod_stan)
                dev.off()
                
                # compile model results
                res_stan_full <- extract(mod_stan)
                summary_beta_stan <- apply(res_stan_full$beta,2,quantile,quantiles)
                
                # summaries of the model parameter estimates
                if (cause_re_model) {
                    summary_sigma_re_stan <- apply(res_stan_full$sigma,2,quantile,quantiles)
                    summary_Omega_stan <- apply(res_stan_full$Omega,c(2,3),quantile,quantiles)
                    cilOmega <- as.vector(summary_Omega_stan[1,,][lower.tri(Omega)])
                    cimOmega <- as.vector(summary_Omega_stan[2,,][lower.tri(Omega)])
                    ciuOmega <- as.vector(summary_Omega_stan[3,,][lower.tri(Omega)])
                    summary_rho_stan <- rbind(cilOmega,cimOmega,ciuOmega)
                    summary_b_stan <- apply(res_stan_full$b,2,quantile,quantiles)
                    dat$b_pred <- summary_b_stan[2,]
                }
                if (time_rw_model) {
                    summary_sigma_rw_stan <- quantile(res_stan_full$sigma_gamma,quantiles)
                    summary_gamma_stan <- apply(res_stan_full$gamma_vec,2,quantile,quantiles)
                    dat$gamma_pred <- summary_gamma_stan[2,]
                }
                if (no_corr_model) {
                    summary_sigma_re_stan <- quantile(res_stan_full$sigma,quantiles)
                    summary_sigma_rw_stan <- quantile(res_stan_full$sigma_gamma,quantiles)
                    summary_b_stan <- apply(res_stan_full$b,2,quantile,quantiles)
                    summary_gamma_stan <- apply(res_stan_full$gamma_vec,2,quantile,quantiles)
                    dat$b_pred <- summary_b_stan[2,]
                    dat$gamma_pred <- summary_gamma_stan[2,]
                }
                
                summary_alpha_stan <- apply(res_stan_full$alpha,2,quantile,quantiles)
                summary_lambda_stan <- apply(res_stan_full$lambda,2,quantile,quantiles)
                
                ## merge model component predictions onto data
                dat$alpha_pred <- summary_alpha_stan[2,]
                dat$lambda_pred <- summary_lambda_stan[2,]
                
                ## plot results
                setwd(savedir)
                pdf(paste0(savedir, "graphs/sims/explore_correlation/model_",model_number,"_re_compare_truth.pdf"),width=6,height=6)
                
                plot(beta[1],type="n",ylim=range(c(summary_beta_stan[,1],beta[1])),xlim=c(0.5,1.5),
                     xlab=expression(alpha),xaxt="n",ylab="value")
                axis(1,at=1,labels=NA)
                abline(h=beta[1],lty=2,col=alpha("red",0.75))
                points(1,summary_beta_stan[2,1],pch=18,col="dodgerblue",cex=1.5)
                segments(1,summary_beta_stan[1,1],1,summary_beta_stan[3,1],col="dodgerblue")
                legend(1,max(c(summary_beta_stan[,1],beta[1]))+0.2*diff(range(c(summary_beta_stan[,1],beta[1]))),c("truth","Stan w/ 80% CI"),
                       col=c("red","dodgerblue"),pch=c(NA,18),lty=c(2,1),ncol=2,bty = "n",xpd=TRUE,cex=0.7)
                
                plot(rep(beta[-1],ncol(X)),type="n",ylim=range(c(summary_beta_stan[,-1],beta[-1])),xlim=c(1,ncol(X)-1),
                     xlab=expression(beta),xaxt="n",ylab="value")
                axis(1,at=(1:ncol(X))+0.3,labels=NA)
                points(1:(ncol(X)-1),beta[-1],pch=19,col=alpha("red",0.75))
                points((1:(ncol(X)-1))+0.2,c(summary_beta_stan[2,-1]),pch=18,col="dodgerblue",cex=1.5)
                segments((1:(ncol(X)-1))+0.2,summary_beta_stan[1,-1],(1:(ncol(X)-1))+0.2,summary_beta_stan[3,-1],col="dodgerblue")
                legend(1,max(c(summary_beta_stan[,-1],beta[-1]))+0.2*diff(range(c(summary_beta_stan[,-1],beta[-1]))),c("truth","Stan w/ 80% CI"),
                       col=c("red","dodgerblue"),pch=c(NA,18),lty=c(2,1),ncol=2,bty = "n",xpd=TRUE,cex=0.7)
                
                if (cause_re_model) {
                    plot(rep(sigma_re,ncause),type="n",ylim=range(c(summary_sigma_re_stan,sigma_re)),xlim=c(1,ncause+0.5),
                     xlab=expression(sigma),xaxt="n",ylab="value")
                    axis(1,at=1:ncause,labels=NA)
                    abline(h=sigma_re,lty=2,col=alpha("red",0.75))
                    points(1:ncause,c(summary_sigma_re_stan[2,]),pch=18,col="dodgerblue",cex=1.5)
                    segments(1:ncause,summary_sigma_re_stan[1,],1:ncause,summary_sigma_re_stan[3,],col="dodgerblue")
                    legend(1,max(c(summary_sigma_re_stan,sigma_re))+0.15*diff(range(c(summary_sigma_re_stan,sigma_re))),c("truth","Stan w/ 80% CI"),
                           col=c("red","dodgerblue"),pch=c(NA,18),lty=c(2,1),ncol=3,bty = "n",xpd=TRUE,cex=0.7)
                }
                
                if (no_corr_model) {
                    plot(rep(sigma_re,ncause),type="n",ylim=range(c(summary_sigma_re_stan,sigma_re)),xlim=c(1,ncause+0.5),
                     xlab=expression(sigma),xaxt="n",ylab="value")
                    axis(1,at=1:ncause,labels=NA)
                    abline(h=sigma_re,lty=2,col=alpha("red",0.75))
                    points(1,c(summary_sigma_re_stan[2]),pch=18,col="dodgerblue",cex=1.5)
                    segments(1,summary_sigma_re_stan[1],1,summary_sigma_re_stan[3],col="dodgerblue")
                    legend(1,max(c(summary_sigma_re_stan,sigma_re))+0.15*diff(range(c(summary_sigma_re_stan,sigma_re))),c("truth","Stan w/ 80% CI"),
                           col=c("red","dodgerblue"),pch=c(NA,18),lty=c(2,1),ncol=3,bty = "n",xpd=TRUE,cex=0.7)
                }
                
                if (time_rw_model | no_corr_model) {
                     plot(sigma_rw,type="n",ylim=range(c(summary_sigma_rw_stan,sigma_rw)),xlim=c(0.5,1.5),
                     xlab=expression(sigma[gamma]),xaxt="n",ylab="value")
                    axis(1,at=1,labels=NA)
                    abline(h=sigma_rw,lty=2,col=alpha("red",0.75))
                    points(1,c(summary_sigma_rw_stan[2]),pch=18,col="dodgerblue",cex=1.5)
                    segments(1,summary_sigma_rw_stan[1],1,summary_sigma_rw_stan[3],col="dodgerblue")
                    legend(1,max(c(summary_sigma_rw_stan,sigma_rw))+0.15*diff(range(c(summary_sigma_rw_stan,sigma_rw))),c("truth","Stan w/ 80% CI"),
                           col=c("red","dodgerblue"),pch=c(NA,18),lty=c(2,1),ncol=3,bty = "n",xpd=TRUE,cex=0.7)
                }
               
                if (cause_re_model) {
                    if (ncause == 2) {
                        plot(rho,type="n",ylim=range(c(summary_rho_stan,rho)),xlim=c(1,length(true_rhos)+0.5),
                             xlab=expression(rho),xaxt="n",ylab="value")
                        axis(1,at=c(1:(length(true_rhos)+0.3)),labels=NA)
                        points(1:length(true_rhos)+0.15,true_rhos,col="red",pch=19)
                        points(1:(length(true_rhos)),c(summary_rho_stan[2]),pch=18,col="dodgerblue",cex=1.5)
                        segments(1:(length(true_rhos)),summary_rho_stan[1],1:(length(true_rhos)),summary_rho_stan[3],col="dodgerblue")
                        legend(1,max(c(summary_rho_stan,rho))+0.15*diff(range(c(summary_rho_stan,rho))),c("truth","Stan w/ 80% CI"),
                               col=c("red","dodgerblue"),pch=c(NA,18),lty=c(2,1),ncol=3,bty = "n",xpd=TRUE,cex=0.7)
                    } else {
                        rho_order <- order(true_rhos)
                        plot(true_rhos,type="n",ylim=range(c(summary_rho_stan,true_rhos)),xlim=c(1,length(true_rhos)+0.5),
                             xlab=expression(rho),xaxt="n",ylab="value")
                        axis(1,at=c(1:(length(true_rhos)+0.3)),labels=NA)
                        points(1:length(true_rhos)+0.15,true_rhos[rho_order],col="red",pch=19)
                        points(1:(length(true_rhos)),c(summary_rho_stan[2,rho_order]),pch=18,col="dodgerblue",cex=1.5)
                        segments(1:(length(true_rhos)),summary_rho_stan[1,rho_order],1:(length(true_rhos)),summary_rho_stan[3,rho_order],col="dodgerblue")
                        legend(1,max(c(summary_rho_stan,rho))+0.15*diff(range(c(summary_rho_stan,rho))),c("truth","Stan w/ 80% CI"),
                               col=c("red","dodgerblue"),pch=c(NA,18),lty=c(2,1),ncol=3,bty = "n",xpd=TRUE,cex=0.7)
                    }
                }
                
                dev.off()
                
                ## plot a visualization of the components of the model
                pdf(paste0(savedir, "graphs/sims/explore_correlation/model_",model_number,"_example_component_plot.pdf"),height=8,width=10)
                
                ## example plot with just one
                exindex <- dat$reg==2 & dat$cause==2 & dat$ageg==2
                
                par(mfcol=c(1,1), mar=c(5.1, 4.1, 2.1, 10.6), xpd=NA)
                
                ylim <- range(dat[exindex,c("logmx","lambda","lambda_pred","alpha","alpha_pred")])
                plot(logmx~year,data=dat[exindex,],col=alpha("purple",0.75),pch=19,
                     ylim=ylim,
                     main="",ylab="log(mortality)",xlab="year")
                lines(lambda~year,data=dat[exindex,],col=alpha("tomato4",0.75))
                segments(0.5,mean(dat$alpha[exindex]),nyear+2,
                         mean(dat$alpha[exindex]),col="tomato1",lty=3)
                if (cause_re_model | no_corr_model) {
                    segments(0.5,mean(dat$b[exindex])+mean(dat$alpha[exindex]),nyear+2,
                             mean(dat$b[exindex])+mean(dat$alpha[exindex]),col="tomato2",lty=2)
                    arrows(nyear+2,mean(dat$alpha[exindex]),nyear+2,mean(dat$b[exindex])+mean(dat$alpha[exindex]),
                           col="tomato2",length=.2,lty=2)
                }
                
                lines(lambda_pred~year,data=dat[exindex,],col=alpha("dodgerblue4",0.75))
                segments(0.5,mean(dat$alpha_pred[exindex]),nyear+3,
                         mean(dat$alpha_pred[exindex]),col="dodgerblue1",lty=3)
                if (cause_re_model | no_corr_model) {
                    segments(0.5,mean(dat$b_pred[exindex])+mean(dat$alpha_pred[exindex]),nyear+3,
                             mean(dat$b_pred[exindex])+mean(dat$alpha_pred[exindex]),col="dodgerblue2",lty=2)
                    arrows(nyear+3,mean(dat$alpha_pred[exindex]),nyear+3,mean(dat$b_pred[exindex])+mean(dat$alpha_pred[exindex]),
                           col="dodgerblue2",length=.2,lty=2)
                }
                
                legend("topright",c("truth","model est","sim data","FE","FE + corr RE","FE + corr RE + RW2"),inset=c(-0.3,0),
                       pch=c(NA,NA,19,NA,NA,NA),lty=c(NA,NA,NA,3,2,1),fill=c("tomato","dodgerblue",NA,NA,NA,NA),
                       col=c(NA,NA,"purple","black","black","black"),bty="n",border=c(NA,NA,NA,NA,NA,NA),cex=1)
                
                dev.off()
                
                plotrows <- 14
                plotcols <- 8
                
                pdf(paste0(savedir, "graphs/sims/explore_correlation/model_",model_number,"_all_component_plots.pdf"),height=plotrows*2,width=plotcols*2)
                
                par(mfcol=c(14,8), mar=c(2.1, 1.1, 1.1, 4.1), xpd=NA,oma=c(4,4,1,1))
                counter <- 0
                for (exreg in 1:nreg) {
                    for (exageg in 1:nageg) {
                        for (excause in 1:ncause) {
                            counter <- counter+1
                            exindex <- dat$reg==exreg & dat$cause==excause & dat$ageg==exageg
                            
                            ylim <- range(dat[exindex,c("logmx","lambda","lambda_pred","alpha","alpha_pred")])
                            plot(logmx~year,data=dat[exindex,],col=alpha("purple",0.75),pch=19,
                                 ylim=ylim,
                                 main="",ylab="",xlab="")
                            if (counter==1) text(-6,mean(ylim),"log(mortality)",srt=90)
                            if (counter==plotrows) text((nyear+2)/2,min(ylim)-diff(ylim)*0.3,"year")
                            lines(lambda~year,data=dat[exindex,],col=alpha("tomato4",0.75))
                            segments(0.5,mean(dat$alpha[exindex]),nyear+2,
                                     mean(dat$alpha[exindex]),col="tomato1",lty=3)
                            if (cause_re_model | no_corr_model) {
                                segments(0.5,mean(dat$b[exindex])+mean(dat$alpha[exindex]),nyear+2,
                                         mean(dat$b[exindex])+mean(dat$alpha[exindex]),col="tomato2",lty=2)
                                arrows(nyear+2,mean(dat$alpha[exindex]),nyear+2,mean(dat$b[exindex])+mean(dat$alpha[exindex]),
                                       col="tomato2",length=.05,lty=2)
                            }
                            lines(lambda_pred~year,data=dat[exindex,],col=alpha("dodgerblue4",0.75))
                            segments(0.5,mean(dat$alpha_pred[exindex]),nyear+3,
                                     mean(dat$alpha_pred[exindex]),col="dodgerblue1",lty=3)
                            if (cause_re_model | no_corr_model) {
                                segments(0.5,mean(dat$b_pred[exindex])+mean(dat$alpha_pred[exindex]),nyear+3,
                                         mean(dat$b_pred[exindex])+mean(dat$alpha_pred[exindex]),col="dodgerblue2",lty=2)
                                arrows(nyear+3,mean(dat$alpha_pred[exindex]),nyear+3,mean(dat$b_pred[exindex])+mean(dat$alpha_pred[exindex]),
                                       col="dodgerblue2",length=.05,lty=2)
                            }
                        }
                    }
                }
                plot.new()
                legend("center",c("truth","model est","sim data","FE","FE + corr RE","FE + corr RE + RW2"),inset=c(-0.3,0),
                       pch=c(NA,NA,19,NA,NA,NA),lty=c(NA,NA,NA,3,2,1),fill=c("tomato","dodgerblue",NA,NA,NA,NA),
                       col=c(NA,NA,"purple","black","black","black"),bty="n",border=c(NA,NA,NA,NA,NA,NA),cex=1.5)
                dev.off()
                
                ## posterior dist of parameters
                if (root != "/home/students/aeschuma/" & ncause == 8 & nreg==6 & nageg==6) {
                    posterior <- as.matrix(mod_stan)
                    
                    plot_title_posterior <- ggtitle("Posterior distributions",
                                                    "with medians and 80% intervals")
                    
                    pdf(paste0(savedir, "graphs/sims/explore_correlation/model_",model_number,"_posteriordists.pdf"), 
                        width=10,height=16)
                    
                    # betas
                    print(mcmc_areas(posterior,
                               pars = c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]","beta[6]","beta[7]","beta[8]",
                                        "beta[9]","beta[10]","beta[11]","beta[12]","beta[13]","beta[14]","beta[15]","beta[16]","beta[17]","beta[18]"),
                               prob = 0.8) + plot_title_posterior)
                    
                    # cause RE correlations
                    if (cause_re_model) {
                        print(mcmc_areas(posterior,
                                   pars = c("Omega[2,1]","Omega[3,1]","Omega[4,1]","Omega[5,1]","Omega[6,1]","Omega[7,1]","Omega[8,1]",
                                            "Omega[3,2]","Omega[4,2]","Omega[5,2]","Omega[6,2]","Omega[7,2]","Omega[8,2]",
                                            "Omega[4,3]","Omega[5,3]","Omega[6,3]","Omega[7,3]","Omega[8,3]",
                                            "Omega[5,4]","Omega[6,4]","Omega[7,4]","Omega[8,4]",
                                            "Omega[6,5]","Omega[7,5]","Omega[8,5]",
                                            "Omega[7,6]","Omega[8,6]",
                                            "Omega[8,7]"),
                                   prob = 0.8) + plot_title_posterior)
                    }
                    # time RW variance
                    # cause RE variances
                    if (cause_re_model) {
                        print(mcmc_areas(posterior,
                                   pars = c("sigma[1]","sigma[2]","sigma[3]","sigma[4]","sigma[5]","sigma[6]","sigma[7]","sigma[8]"),
                                   prob = 0.8) + plot_title_posterior)
                    }
                    if (no_corr_model) {
                        print(mcmc_areas(posterior,
                                   pars = c("sigma"),
                                   prob = 0.8) + plot_title_posterior)
                    }
                    # for general exam presentation
                    # mcmc_areas(posterior,
                    #            pars = c("Omega[2,1]","Omega[6,7]","Omega[4,1]"),
                    #            prob = 0.8) + plot_title_posterior
                    if (time_rw_model | no_corr_model) {
                        print(mcmc_areas(posterior,
                                   pars = c("sigma_gamma"),
                                   prob = 0.8) + plot_title_posterior)
                    }
                    dev.off()
                }
            }
        }
    }
}

