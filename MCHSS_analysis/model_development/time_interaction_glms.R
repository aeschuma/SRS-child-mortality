## Austin Schumacher
## 12/17/2019
## Fit a bunch of GLMs to MCHSS data and explore residuals

# preamble ####
rm(list=ls())

library(tidyverse);library(data.table);library(gtools);
library(RColorBrewer);library(viridis);library(lattice);library(scales);
library(MASS);library(lme4);library(glmmTMB);

# directory for results (set this yourself if you need to change it!)
savedir <- "../../../../Dropbox/SRS-child-mortality-output/"

# create folders to store results if necessary
if (!file.exists(paste0(savedir, "graphs"))) {
    dir.create(paste0(savedir, "graphs"))
}

if (!file.exists(paste0(savedir, "graphs/model_testing"))) {
    dir.create(paste0(savedir, "graphs/model_testing"))
}

# code running options ####
fe_models <- TRUE
time_models <- TRUE

# plotting functions ####

# plotting function for time residuals
plot_time_resids <- function(mymodel,mydata,
                             agecols,reg2cols,causecols,
                             agepch,reg2pch,causepch,
                             agelty,reg2lty,causelty) {
    
    ## testing
    # mymodel <- glm(deaths ~ offset(logpy) +
    #                    factor(reg2) * factor(agegp) +
    #                    factor(reg2) * factor(cause) +
    #                    factor(agegp) * factor(cause) +
    #                    (year + year2)*factor(agereg2cause),
    #                data=chn_all,family=quasipoisson)
    # mydata <- chn_all
    # agecols <- myagecols
    # reg2cols <- myreg2cols
    # causecols <- mycausecols
    # agepch <- myagepch
    # reg2pch <- myreg2pch
    # causepch <- mycausepch
    # agelty <- myagelty
    # reg2lty <- myreg2lty
    # causelty <- mycauselty
    
    ## Begin function code
    resids <- resid(mymodel)
    
    plot.new()
    title(mymodel$call,cex.main=0.5)
    
    # plot axis settings
    font.lab <- 2
    cex.main <- 1.5
    cex.lab <- 1.7
    cex.axis <- 1.5
    cex.legend.big <- 1.2
    cex.legend.small <- 1
    
    # plots by age, region, and cause separately
    plot(resids~jitter(mydata$year),col=alpha(agecols,0.3),pch=19,main="Residuals as a function of year \n colored by age",
         xlab="Year (0 refers to 2005.5)",ylab="Residuals",
         font.lab = font.lab, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    legend("top",legend=unique(mydata$agegp_name),col=unique(agecols),pch=19,
           ncol=length(unique(mydata$agegp_name))/2,lty=1,text.font = 2,cex=cex.legend.big)
    for (i in 1:length(unique(mydata$agegp_name))) {
        lines(lowess(resids[mydata$agegp_name==unique(mydata$agegp_name)[i]]~
                         mydata$year[mydata$agegp_name==unique(mydata$agegp_name)[i]]),
              col=agecols[mydata$agegp_name==unique(mydata$agegp_name)[i]][1],lwd=2)
        
    }
    plot(resids~jitter(mydata$year),col=alpha(reg2cols,0.3),pch=19,main="Residuals as a function of year \n colored by region",
         xlab="Year (0 refers to 2005.5)",ylab="Residuals",
         font.lab = font.lab, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    legend("top",legend=unique(mydata$reg2),col=unique(reg2cols),pch=19,
           ncol=length(unique(mydata$reg2))/2,lty=1,text.font = 2,cex=cex.legend.big)
    for (i in 1:length(unique(mydata$reg2))) {
        lines(lowess(resids[mydata$reg2==unique(mydata$reg2)[i]]~
                         mydata$year[mydata$reg2==unique(mydata$reg2)[i]]),
              col=reg2cols[mydata$reg2==unique(mydata$reg2)[i]][1],lwd=2)
        
    }
    plot(resids~jitter(mydata$year),col=alpha(causecols,0.3),pch=19,main="Residuals as a function of year \n colored by cause",
         xlab="Year (0 refers to 2005.5)",ylab="Residuals",
         font.lab = font.lab, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    legend("top",legend=unique(mydata$cause_name),col=unique(causecols),
           pch=19,ncol=length(unique(mydata$cause_name))/4,lty=1,text.font = 2,cex=cex.legend.big)
    for (i in 1:length(unique(mydata$cause_name))) {
        lines(lowess(resids[mydata$cause_name==unique(mydata$cause_name)[i]]~
                         mydata$year[mydata$cause_name==unique(mydata$cause_name)[i]]),
              col=causecols[mydata$cause_name==unique(mydata$cause_name)[i]][1],lwd=2)
        
    }
    
    # plots colored and symboled by age, region, and cause interactions
    plot(resids~jitter(mydata$year),col=alpha(agecols,0.3),pch=reg2pch,
         main="Residuals as a function of year \n colored by age, symboled by region",
         xlab="Year (0 refers to 2005.5)",ylab="Residuals",
         font.lab = font.lab, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    legend("top",legend=c(unique(as.character(mydata$agegp_name)),unique(as.character(mydata$reg2))),
           col=c(unique(agecols),rep("black",length(unique(mydata$reg2)))),
           pch=c(rep(NA,length(unique(mydata$agegp_name))),unique(reg2pch)),
           lty=c(rep(NA,length(unique(mydata$agegp_name))),unique(reg2lty)),
           fill=c(unique(agecols),rep(NA,length(unique(mydata$reg2)))),
           ncol=4,text.font = 2,cex=cex.legend.big,
           border=NA)
    for (i in 1:length(unique(mydata$agereg2))) {
        lines(lowess(resids[mydata$agereg2==unique(mydata$agereg2)[i]]~
                         mydata$year[mydata$agereg2==unique(mydata$agereg2)[i]]),
              col=agecols[mydata$agereg2==unique(mydata$agereg2)[i]][1],
              lty=reg2lty[mydata$agereg2==unique(mydata$agereg2)[i]][1],
              lwd=2)
    }
    plot(resids~jitter(mydata$year),col=alpha(reg2cols,0.3),pch=agepch,
         main="Residuals as a function of year \n colored by region, symboled by age",
         xlab="Year (0 refers to 2005.5)",ylab="Residuals",
         font.lab = font.lab, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    legend("top",legend=c(unique(as.character(mydata$reg2)),unique(as.character(mydata$agegp_name))),
           col=c(unique(reg2cols),rep("black",length(unique(mydata$agegp_name)))),
           pch=c(rep(NA,length(unique(mydata$reg2))),unique(agepch)),
           lty=c(rep(NA,length(unique(mydata$reg2))),unique(agelty)),
           fill=c(unique(reg2cols),rep(NA,length(unique(mydata$agegp_name)))),
           ncol=4,text.font = 2,cex=cex.legend.big,
           border=NA)
    for (i in 1:length(unique(mydata$agereg2))) {
        lines(lowess(resids[mydata$agereg2==unique(mydata$agereg2)[i]]~
                         mydata$year[mydata$agereg2==unique(mydata$agereg2)[i]]),
              col=reg2cols[mydata$agereg2==unique(mydata$agereg2)[i]][1],
              lty=agelty[mydata$agereg2==unique(mydata$agereg2)[i]][1],
              lwd=2)
    }
    
    plot(resids~jitter(mydata$year),col=alpha(agecols,0.3),pch=causepch,
         main="Residuals as a function of year \n colored by age, symboled by cause",
         xlab="Year (0 refers to 2005.5)",ylab="Residuals",
         font.lab = font.lab, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    legend("top",legend=c(unique(as.character(mydata$agegp_name)),unique(as.character(mydata$cause_name))),
           col=c(unique(agecols),rep("black",length(unique(mydata$cause_name)))),
           pch=c(rep(NA,length(unique(mydata$agegp_name))),unique(causepch)),
           lty=c(rep(NA,length(unique(mydata$agegp_name))),unique(causelty)),
           fill=c(unique(agecols),rep(NA,length(unique(mydata$cause_name)))),
           ncol=2,cex=cex.legend.small,text.font = 2,
           border=NA)
    for (i in 1:length(unique(mydata$agecause))) {
        lines(lowess(resids[mydata$agecause==unique(mydata$agecause)[i]]~
                         mydata$year[mydata$agecause==unique(mydata$agecause)[i]]),
              col=agecols[mydata$agecause==unique(mydata$agecause)[i]][1],
              lty=causelty[mydata$agecause==unique(mydata$agecause)[i]][1],
              lwd=2)
    }
    plot(resids~jitter(mydata$year),col=alpha(causecols,0.3),pch=agepch,
         main="Residuals as a function of year \n colored by cause, symboled by age",
         xlab="Year (0 refers to 2005.5)",ylab="Residuals",
         font.lab = font.lab, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    legend("top",legend=c(unique(as.character(mydata$cause_name)),unique(as.character(mydata$agegp_name))),
           col=c(unique(causecols),rep("black",length(unique(mydata$agegp_name)))),
           pch=c(rep(NA,length(unique(mydata$cause_name))),unique(agepch)),
           lty=c(rep(NA,length(unique(mydata$cause_name))),unique(agelty)),
           fill=c(unique(causecols),rep(NA,length(unique(mydata$agegp_name)))),
           ncol=2,cex=cex.legend.small,text.font = 2,
           border=NA)
    for (i in 1:length(unique(mydata$agecause))) {
        lines(lowess(resids[mydata$agecause==unique(mydata$agecause)[i]]~
                         mydata$year[mydata$agecause==unique(mydata$agecause)[i]]),
              col=causecols[mydata$agecause==unique(mydata$agecause)[i]][1],
              lty=agelty[mydata$agecause==unique(mydata$agecause)[i]][1],
              lwd=2)
    }
    
    plot(resids~jitter(mydata$year),col=alpha(reg2cols,0.3),pch=causepch,
         main="Residuals as a function of year \n colored by region, symboled by cause",
         xlab="Year (0 refers to 2005.5)",ylab="Residuals",
         font.lab = font.lab, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    legend("top",legend=c(unique(as.character(mydata$reg2)),unique(as.character(mydata$cause_name))),
           col=c(unique(reg2cols),rep("black",length(unique(mydata$cause_name)))),
           pch=c(rep(NA,length(unique(mydata$reg2))),unique(causepch)),
           lty=c(rep(NA,length(unique(mydata$reg2))),unique(causelty)),
           fill=c(unique(reg2cols),rep(NA,length(unique(mydata$cause_name)))),
           ncol=2,cex=cex.legend.small,text.font = 2,
           border=NA)
    for (i in 1:length(unique(mydata$reg2cause))) {
        lines(lowess(resids[mydata$reg2cause==unique(mydata$reg2cause)[i]]~
                         mydata$year[mydata$reg2cause==unique(mydata$reg2cause)[i]]),
              col=reg2cols[mydata$reg2cause==unique(mydata$reg2cause)[i]][1],
              lty=causelty[mydata$reg2cause==unique(mydata$reg2cause)[i]][1],
              lwd=2)
    }
    plot(resids~jitter(mydata$year),col=alpha(causecols,0.3),pch=reg2pch,
         main="Residuals as a function of year \n colored by cause, symboled by region",
         xlab="Year (0 refers to 2005.5)",ylab="Residuals",
         font.lab = font.lab, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    legend("top",legend=c(unique(as.character(mydata$cause_name)),unique(as.character(mydata$reg2))),
           col=c(unique(causecols),rep("black",length(unique(mydata$reg2)))),
           pch=c(rep(NA,length(unique(mydata$cause_name))),unique(reg2pch)),
           lty=c(rep(NA,length(unique(mydata$cause_name))),unique(reg2lty)),
           fill=c(unique(causecols),rep(NA,length(unique(mydata$reg2)))),
           ncol=2,cex=cex.legend.small,text.font = 2,
           border=NA)
    for (i in 1:length(unique(mydata$reg2cause))) {
        lines(lowess(resids[mydata$reg2cause==unique(mydata$reg2cause)[i]]~
                         mydata$year[mydata$reg2cause==unique(mydata$reg2cause)[i]]),
              col=causecols[mydata$reg2cause==unique(mydata$reg2cause)[i]][1],
              lty=reg2lty[mydata$reg2cause==unique(mydata$reg2cause)[i]][1],
              lwd=2)
    }
}

# load and format data ####
chn_all <- as.data.table(readRDS("../../data/mchss_test_data.RDS"))

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

## interaction variables
chn_all$agereg2 <- interaction(chn_all$agegp_name,chn_all$reg2)
chn_all$agecause <- interaction(chn_all$agegp_name,chn_all$cause_name)
chn_all$reg2cause <- interaction(chn_all$reg2,chn_all$cause_name)
chn_all$agereg2cause <- interaction(chn_all$agegp_name,chn_all$reg2,chn_all$cause_name)
chn_all$rw_index <- as.numeric(interaction(chn_all$ageRW,chn_all$reg2,chn_all$causeRW))

# Modeling ####

## FE interaction modeling
if (fe_models) {
    mod1 <- glm(deaths ~ offset(logpy) + factor(agegp) + factor(reg2) + factor(cause) + year,data=chn_all,family=poisson)
    mod2 <- glm(deaths ~ offset(logpy) + factor(agegp) * factor(reg2) * factor(cause) + year,data=chn_all,family=poisson)
    mod3 <- glm(deaths ~ offset(logpy) + factor(reg2) + factor(agegp) * factor(cause) + year,data=chn_all,family=poisson)
    mod4 <- glm(deaths ~ offset(logpy) + factor(reg2) * factor(agegp) + factor(cause) + year,data=chn_all,family=poisson)
    mod5 <- glm(deaths ~ offset(logpy) + factor(agegp) + factor(reg2) * factor(cause) + year,data=chn_all,family=poisson)
    mod6 <- glm(deaths ~ offset(logpy) + factor(reg2) * factor(agegp) + 
                    factor(reg2) * factor(cause) + 
                    factor(agegp) * factor(cause) + 
                    year,
                data=chn_all,family=poisson)
    AIC(mod1,mod2,mod3,mod4,mod5,mod6)
}

## Time modeling
if (time_models) {
    
    # plotting colors
    myagecols <- brewer.pal(length(unique(chn_all$agegp_name)),"Set1")[as.numeric(as.factor(chn_all$agegp_name))]
    myreg2cols <- brewer.pal(length(unique(chn_all$reg2)),"Set1")[as.numeric(as.factor(chn_all$reg2))]
    mycausecols <- brewer.pal(length(unique(chn_all$cause_name)),"Set1")[as.numeric(as.factor(chn_all$cause_name))]
    
    # plotting pch
    myagepch <- c(15,17,18,19,7,8)[as.numeric(as.factor(chn_all$agegp_name))]
    myreg2pch <- c(15,17,18,19,7,8)[as.numeric(as.factor(chn_all$reg2))]
    mycausepch <- c(15,17,18,19,7,8,11,3)[as.numeric(as.factor(chn_all$cause_name))]
    
    # plotting lty
    myagelty <- c(1,2,3,4,5,6)[as.numeric(as.factor(chn_all$agegp_name))]
    myreg2lty <- c(1,2,3,4,5,6)[as.numeric(as.factor(chn_all$reg2))]
    mycauselty <- c("solid","dashed","dotted","dotdash","longdash","twodash","1F","12345678")[as.numeric(as.factor(chn_all$cause_name))]
    
    # overall linear and quadratic time effect
    tmod1 <- glm(deaths ~ offset(logpy) + 
                     factor(reg2) * factor(agegp) + 
                     factor(reg2) * factor(cause) + 
                     factor(agegp) * factor(cause) + 
                     year + year2,
                 data=chn_all,family=quasipoisson)
    
    # plot residuals
    mywidth <- 8
    myheight <- 9
    
    pdf(paste0(savedir,"graphs/model_testing/time_diag_resids_overall_time.pdf"),
        width=mywidth,height=myheight)
    
    plot_time_resids(tmod1,chn_all,
                     myagecols,myreg2cols,mycausecols,
                     myagepch,myreg2pch,mycausepch,
                     myagelty,myreg2lty,mycauselty)
    
    dev.off()
    
    # linear and quadratic time effects by region
    tmod2 <- glm(deaths ~ offset(logpy) + 
                     factor(reg2) * factor(agegp) + 
                     factor(reg2) * factor(cause) + 
                     factor(agegp) * factor(cause) + 
                     (year + year2)*factor(reg2),
                 data=chn_all,family=quasipoisson)
    
    # plot residuals
    pdf(paste0(savedir,"graphs/model_testing/time_diag_resids_time_by_region.pdf"),
        width=mywidth,height=myheight)
    
    plot_time_resids(tmod2,chn_all,
                     myagecols,myreg2cols,mycausecols,
                     myagepch,myreg2pch,mycausepch,
                     myagelty,myreg2lty,mycauselty)
    
    dev.off()
    
    # linear and quadratic time effects by region
    tmod3 <- glm(deaths ~ offset(logpy) + 
                     factor(reg2) * factor(agegp) + 
                     factor(reg2) * factor(cause) + 
                     factor(agegp) * factor(cause) + 
                     (year + year2)*factor(agegp),
                 data=chn_all,family=quasipoisson)
    
    # plot residuals
    pdf(paste0(savedir,"graphs/model_testing/time_diag_resids_time_by_age.pdf"),
        width=mywidth,height=myheight)
    
    plot_time_resids(tmod3,chn_all,
                     myagecols,myreg2cols,mycausecols,
                     myagepch,myreg2pch,mycausepch,
                     myagelty,myreg2lty,mycauselty)
    
    dev.off()
    
    # linear and quadratic time effects by region
    tmod4 <- glm(deaths ~ offset(logpy) + 
                     factor(reg2) * factor(agegp) + 
                     factor(reg2) * factor(cause) + 
                     factor(agegp) * factor(cause) + 
                     (year + year2)*factor(cause),
                 data=chn_all,family=quasipoisson)
    
    # plot residuals
    pdf(paste0(savedir,"graphs/model_testing/time_diag_resids_time_by_cause.pdf"),
        width=mywidth,height=myheight)
    
    plot_time_resids(tmod4,chn_all,
                     myagecols,myreg2cols,mycausecols,
                     myagepch,myreg2pch,mycausepch,
                     myagelty,myreg2lty,mycauselty)
    
    dev.off()
    
    # linear and quadratic time effects by age, region
    tmod5 <- glm(deaths ~ offset(logpy) + 
                     factor(reg2) * factor(agegp) + 
                     factor(reg2) * factor(cause) + 
                     factor(agegp) * factor(cause) + 
                     (year + year2)*factor(agereg2),
                 data=chn_all,family=quasipoisson)
    
    # plot residuals
    pdf(paste0(savedir,"graphs/model_testing/time_diag_resids_time_by_agereg2.pdf"),
        width=mywidth,height=myheight)
    
    plot_time_resids(tmod5,chn_all,
                     myagecols,myreg2cols,mycausecols,
                     myagepch,myreg2pch,mycausepch,
                     myagelty,myreg2lty,mycauselty)
    
    dev.off()
    
    # linear and quadratic time effects by age, cause
    tmod6 <- glm(deaths ~ offset(logpy) + 
                     factor(reg2) * factor(agegp) + 
                     factor(reg2) * factor(cause) + 
                     factor(agegp) * factor(cause) + 
                     (year + year2)*factor(agecause),
                 data=chn_all,family=quasipoisson)
    
    # plot residuals
    pdf(paste0(savedir,"graphs/model_testing/time_diag_resids_time_by_agecause.pdf"),
        width=mywidth,height=myheight)
    
    plot_time_resids(tmod6,chn_all,
                     myagecols,myreg2cols,mycausecols,
                     myagepch,myreg2pch,mycausepch,
                     myagelty,myreg2lty,mycauselty)
    
    dev.off()
    
    # linear and quadratic time effects by age, region
    tmod7 <- glm(deaths ~ offset(logpy) + 
                     factor(reg2) * factor(agegp) + 
                     factor(reg2) * factor(cause) + 
                     factor(agegp) * factor(cause) + 
                     (year + year2)*factor(reg2cause),
                 data=chn_all,family=quasipoisson)
    
    # plot residuals
    pdf(paste0(savedir,"graphs/model_testing/time_diag_resids_time_by_reg2cause.pdf"),
        width=mywidth,height=myheight)
    
    plot_time_resids(tmod7,chn_all,
                     myagecols,myreg2cols,mycausecols,
                     myagepch,myreg2pch,mycausepch,
                     myagelty,myreg2lty,mycauselty)
    
    dev.off()
    
    # for fake data, need to delete groups with all 0 deaths
    summary_dat <- chn_all %>% group_by(reg2, agegp, cause) %>% 
        summarise(min_deaths = min(deaths),
                  med_deaths = median(deaths),
                  max_deaths = max(deaths),
                  all.zero = sum(deaths) == 0)
    
    chn_all_tmp <- chn_all %>%
        merge(summary_dat)  %>%
        filter(!all.zero)
    
    # linear and quadratic time effects by age, region
    tmod8 <- glm(deaths ~ offset(logpy) + 
                     factor(reg2) * factor(agegp) + 
                     factor(reg2) * factor(cause) + 
                     factor(agegp) * factor(cause) + 
                     (year + year2)*factor(agereg2cause),
                 data=chn_all_tmp,family=quasipoisson)
    
    # plot residuals
    pdf(paste0(savedir,"graphs/model_testing/time_diag_resids_time_by_agereg2cause.pdf"),
        width=mywidth,height=myheight)
    
    plot_time_resids(tmod8,chn_all_tmp,
                     myagecols,myreg2cols,mycausecols,
                     myagepch,myreg2pch,mycausepch,
                     myagelty,myreg2lty,mycauselty)
    
    dev.off()
}
