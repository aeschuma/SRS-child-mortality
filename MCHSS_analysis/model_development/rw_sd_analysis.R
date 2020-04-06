## Austin Schumacher
## 12/19/2019
## Using China data
## Fit a RW model in INLA for each strata we want to do a RW on
## plot the estimated RW SD parameters for comparison

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
library(INLA); library(scales); library(RColorBrewer);

########################
## code running options
########################

#####################
## model options
#####################
quantiles <- c(0.1,0.5,0.9)

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

# load test data
chn_all <- readRDS("../china_data/mchss_test_data.RDS")

## format data
chn_all$logpy <- log(chn_all$exposure)
chn_all$year2 <- chn_all$year^2
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

############
## RW modeling (this also plots the resulting predictions w/ the data)
############

## set up storage for predicted RW SD
rw_sds <- chn_all[,c("ageRW","reg2","causeRW","rw_index")]
rw_sds <- rw_sds[!duplicated(rw_sds),]
rw_sds <- rw_sds[order(rw_sds$rw_index),]
rw_sds$strata_name <- paste(rw_sds$ageRW,rw_sds$reg2,rw_sds$causeRW,sep="_")
rw_sds$sd_med <- rw_sds$sd_lower <- rw_sds$sd_upper <- NA

pdf("graphs/rw_analysis/preds.pdf",
    width=8,height=8)

for (i in 1:nrow(rw_sds)) {
    ## testing
    # i <- 3
    
    ## subset data
    tmpdat <- chn_all[chn_all$rw_index==i,]
    
    ## RW model
    quantile_names <- paste0(quantiles,"quant")
    if (length(unique(tmpdat$agereg2cause)) == 1) {
        mod_inla_test <- inla(deaths ~ 1 + 
                                  f(year, model="rw2", constr=TRUE) +
                                  f(factor(obs), model="iid") +
                                  offset(logpy),
                              family="poisson", 
                              data=tmpdat,
                              quantiles = quantiles,
                              control.predictor=list(compute=TRUE))
        ## FE summary
        summary_beta <- unlist(mod_inla_test$summary.fixed[,quantile_names])
        tmpdat$pred_beta <- summary_beta[2]
    } else {
        mod_inla_test <- inla(deaths ~ factor(agereg2cause) + 
                                  f(year, model="rw2", constr=TRUE) +
                                  f(factor(obs), model="iid") +
                                  offset(logpy),
                              family="poisson", 
                              data=tmpdat,
                              quantiles = quantiles,
                              control.predictor=list(compute=TRUE))
        
        ## FE summary
        summary_beta <- mod_inla_test$summary.fixed[,quantile_names]
        X_tmp <- model.matrix(~factor(agereg2cause),data=tmpdat)
        tmpdat$pred_beta <- as.vector(X_tmp %*% summary_beta[,"0.5quant"])
    }
    
    ## IID summary
    summary_tau_iid <- unlist(mod_inla_test$summary.hyperpar[2,quantile_names])
    summary_sd_iid <- rev(summary_tau_iid^(-0.5))
    tmpdat$pred_iid <- mod_inla_test$summary.random$`factor(obs)`[,"0.5quant"]
    
    ## RW summary
    summary_tau_rw <- unlist(mod_inla_test$summary.hyperpar[1,quantile_names])
    summary_sd_rw <- rev(summary_tau_rw^(-0.5))
    pred_rw <- mod_inla_test$summary.random$year[,c("ID","0.5quant")]
    names(pred_rw) <- c("year","pred_rw")
    tmpdat <- merge(tmpdat,pred_rw)
    tmpdat$preds <- tmpdat$pred_beta + tmpdat$pred_iid + tmpdat$pred_rw
    
    ## plotting settings
    tmpdat$strata <- as.numeric(factor(as.character(tmpdat$agereg2cause)))
    plotcols <- brewer.pal(length(unique(tmpdat$agereg2cause)),"Set1")
    plotpch <- c(0,1,2,5,6,7,8,9,10,11)[1:length(unique(tmpdat$agereg2cause))]
    tmpdat$plotpch <- plotpch[tmpdat$strata]
    
    ## plotting adjustments
    smallest_logmx <- min(c(tmpdat$preds,tmpdat$logmx[!is.infinite(tmpdat$logmx)]))
    tmpdat$logmx_plot <- ifelse(is.infinite(tmpdat$logmx),smallest_logmx-1,tmpdat$logmx)
    tmpdat$logmx_pch <- ifelse(is.infinite(tmpdat$logmx),17,19)
    
    plot(logmx~year,data=tmpdat,type="n",
         ylim=range(c(tmpdat$preds,tmpdat$logmx_plot)),
         main=paste0("age group ",rw_sds$ageRW[rw_sds$rw_index==i],
                     "; ",rw_sds$reg2[rw_sds$rw_index==i],
                     "; ",rw_sds$causeRW[rw_sds$rw_index==i]))
    
    for (j in 1:length(unique(tmpdat$agereg2cause))) {
        plotdat <- tmpdat[tmpdat$strata==j,]
        points(logmx_plot~year,data=plotdat,
               col=plotcols[j],
               pch=plotdat$logmx_pch)
        points(preds~year,data=plotdat,
               col=plotcols[j],
               pch=plotpch[j])
        abline(h=plotdat$pred_beta[1],col=alpha(plotcols[j],0.5))
        lines((plotdat$pred_rw + plotdat$pred_beta)~plotdat$year,
              col=alpha(plotcols[j],0.5),lty=2)
    }
    
    ## store predicted RW SD
    rw_sds$sd_med[rw_sds$rw_index==i] <- summary_sd_rw[2]
    rw_sds$sd_lower[rw_sds$rw_index==i] <- summary_sd_rw[1]
    rw_sds$sd_upper[rw_sds$rw_index==i] <- summary_sd_rw[3]
}
dev.off()

##################
## Plot the predicted RW SDs
##################

# cause colors
cause_cols <- brewer.pal(length(unique(rw_sds$causeRW)),"Set2")
rw_sds$cause_col <- cause_cols[factor(rw_sds$causeRW)]

# region colors
reg2_cols <- brewer.pal(length(unique(rw_sds$reg2)),"Set2")
rw_sds$reg2_col <- reg2_cols[factor(rw_sds$reg2)]

# age colors
age_cols <- brewer.pal(length(unique(rw_sds$ageRW)),"Set1")
rw_sds$age_col <- age_cols[factor(rw_sds$ageRW)]

# region pchs
reg2_pchs <- letters[1:length(unique(rw_sds$reg2))]
rw_sds$reg2_pch <- reg2_pchs[factor(rw_sds$reg2)]

# age pchs
rw_sds$age_pch <- as.character(rw_sds$ageRW)

# cause pchs when paired with age
cause_pchs <- letters[1:length(unique(rw_sds$causeRW))]
rw_sds$cause_pch <- cause_pchs[factor(rw_sds$causeRW)]

# cause pchs when paired with region
cause_pchs2 <- as.character(1:length(unique(rw_sds$causeRW)))
rw_sds$cause_pch2 <- cause_pchs2[factor(rw_sds$causeRW)]

pdf("graphs/rw_analysis/rw_sds.pdf",
    width=16,height=9)

# graphical parameters
par(mar = c(5,5,4,2))

# ordered by cause, age, region
rw_sds <- rw_sds[order(rw_sds$causeRW,rw_sds$ageRW),]
rw_sds$rw_index <- 1:nrow(rw_sds)

# plotting locations (offset for points)
loc.adj <- 0.3
rw_sds$loc1 <- rw_sds$rw_index-loc.adj
rw_sds$loc2 <- rw_sds$rw_index+loc.adj

cex <- 1.5
cex.main <- 1.7
cex.lab <- 1.6
cex.axis <- 1.5
cex.legend <- 1.8

plot(rw_sds$sd_med~rw_sds$loc1,
     xlim=c(1,nrow(rw_sds)),
     ylim=range(c(rw_sds$sd_lower,rw_sds$sd_upper+0.03)),
     pch=rw_sds$reg2_pch,
     col=rw_sds$cause_col,
     xlab="predicted RW SD",ylab="RW index",
     main="Posterior median RW SDs \n w/ 80% posterior predictive interval",
     cex=cex,cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis,font.lab=2)
points(rw_sds$sd_med~rw_sds$loc2,
       pch=rw_sds$age_pch,
       col=rw_sds$cause_col,cex=cex)
segments(rw_sds$rw_index,rw_sds$sd_lower,rw_sds$rw_index,rw_sds$sd_upper,
         col=rw_sds$cause_col)
legend(25,0.4,c(unique(as.character(factor(rw_sds$causeRW))),
                    unique(as.character(factor(rw_sds$reg2))),
                    "neonates","older ages"),
       col=c(cause_cols,rep("black",length(reg2_pchs)),rep("black",2)),
       pch=c(rep(NA,length(cause_cols)),reg2_pchs,"1","2"),
       fill=c(cause_cols,rep(NA,length(reg2_pchs)),rep(NA,2)),
       border=NA,ncol=2,cex=cex.legend,bty='n',text.font = 2)

# ordered by region, age, cause
rw_sds <- rw_sds[order(rw_sds$reg2,rw_sds$ageRW),]
rw_sds$rw_index <- 1:nrow(rw_sds)

# plotting locations (offset for points)
rw_sds$loc1 <- rw_sds$rw_index-loc.adj
rw_sds$loc2 <- rw_sds$rw_index+loc.adj

plot(rw_sds$sd_med~rw_sds$loc1,
     xlim=c(1,nrow(rw_sds)),
     ylim=range(c(rw_sds$sd_lower,rw_sds$sd_upper+0.03)),
     pch=rw_sds$cause_pch,
     col=rw_sds$reg2_col,
     xlab="predicted RW SD",ylab="RW index",
     main="Posterior median RW SDs \n w/ 80% posterior predictive interval",
     cex=cex,cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis,font.lab=2)
points(rw_sds$sd_med~rw_sds$loc2,
       pch=rw_sds$age_pch,
       col=rw_sds$reg2_col,cex=cex)
segments(rw_sds$rw_index,rw_sds$sd_lower,rw_sds$rw_index,rw_sds$sd_upper,
         col=rw_sds$reg2_col)
legend("topleft",c(unique(as.character(factor(rw_sds$causeRW))),
               unique(as.character(factor(rw_sds$reg2))),
               "neonates","older ages"),
       col=c(rep("black",length(cause_pchs)),reg2_cols,rep("black",2)),
       pch=c(cause_pchs,rep(NA,length(reg2_cols)),"1","2"),
       fill=c(rep(NA,length(cause_pchs)),reg2_cols,rep(NA,2)),
       border=NA,ncol=2,cex=cex.legend,bty='n',text.font = 2)

# order by age, region, cause
rw_sds <- rw_sds[order(rw_sds$ageRW,rw_sds$reg2),]
rw_sds$rw_index <- 1:nrow(rw_sds)

# plotting locations (offset for points)
rw_sds$loc1 <- rw_sds$rw_index-loc.adj
rw_sds$loc2 <- rw_sds$rw_index+loc.adj

plot(rw_sds$sd_med~rw_sds$loc1,
     xlim=c(1,nrow(rw_sds)),
     ylim=range(c(rw_sds$sd_lower,rw_sds$sd_upper+0.03)),
     pch=rw_sds$reg2_pch,
     col=rw_sds$age_col,
     xlab="predicted RW SD",ylab="RW index",
     main="Posterior median RW SDs \n w/ 80% posterior predictive interval",
     cex=cex,cex.main=cex.main,cex.lab=cex.lab,cex.axis=cex.axis,font.lab=2)
points(rw_sds$sd_med~rw_sds$loc2,
       pch=rw_sds$cause_pch2,
       col=rw_sds$age_col,cex=cex)
segments(rw_sds$rw_index,rw_sds$sd_lower,rw_sds$rw_index,rw_sds$sd_upper,
         col=rw_sds$age_col)
legend(-0.75,0.4,c(unique(as.character(factor(rw_sds$causeRW))),
               unique(as.character(factor(rw_sds$reg2))),
               "neonates","older ages"),
       col=c(rep("black",length(cause_pchs2)),rep("black",length(reg2_pchs)),age_cols[1:2]),
       pch=c(cause_pchs2,reg2_pchs,NA,NA),
       fill=c(rep(NA,length(cause_pchs2)),rep(NA,length(reg2_pchs)),age_cols[1:2]),
       border=NA,ncol=2,cex=cex.legend,bty='n',text.font = 2)

dev.off()
