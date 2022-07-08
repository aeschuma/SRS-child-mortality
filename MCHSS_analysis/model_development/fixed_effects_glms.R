## Austin Schumacher
## 2018 Nov 23
## 1. Fit GLMs on China data to assess the drivers of variability and structure/correlations in residuals
## 2. Fit GLMMs on China data to assess model fit, drivers of variability, and correlation in residuals

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

if (!file.exists(paste0(savedir, "graphs/fe_re_glms"))) {
    dir.create(paste0(savedir, "graphs/fe_re_glms"))
}

# load and format data ####
chn_all <- as.data.table(readRDS("../../data/mchss_test_data.RDS"))

# indices and formatting
causes <- unique(chn_all$cause)
chn_all$logpy <- log(chn_all$exposure)

# All-cause mortality analyses ####

## GLMs ####
dims <- c("year","agegp_name","reg","res")

loops <- list()
for (rr in 1:length(dims)) {
    loops[[rr]] <- combinations(n=length(dims),r=rr,v=dims)
}
mod_glm <- list()
mod_glm_nb <- list()
formulas <- list()
i <- 0
for (r in 1:length(loops)) {
    for (rr in 1:nrow(loops[[r]])) {
        i <- i+1
        formulas[[i]] <- as.formula(paste("deaths ~ offset(logpy) + ",paste(loops[[r]][rr,],collapse=" + "),sep=""))
        mod_glm[[i]] <- glm(formulas[[i]],data=chn_all,family=quasipoisson(link="log"))
        mod_glm_nb[[i]] <- glm.nb(formulas[[i]],data=chn_all)
    }
}

##  AIC analysis
lapply(mod_glm_nb,AIC)
best.glm.nb <- mod_glm_nb[[which(unlist(lapply(mod_glm_nb,AIC))==min(unlist(lapply(mod_glm_nb,AIC))))]]

## plot FEs and residuals
pdf(paste0(savedir,"graphs/fe_re_glms/all_cause_poisson_glm_fe_plots.pdf"),width=16,height=9)
for (i in 1:length(mod_glm)) {
    ncoef <- length(mod_glm[[i]]$coefficients)
    se <- summary(mod_glm[[i]])$coef[,2]
    plot(mod_glm[[i]]$coefficients,
         ylim=range(c(mod_glm[[i]]$coefficients-1.96*se,mod_glm[[i]]$coefficients+1.96*se)),
         xaxt="n",
         main=paste("Poisson GLM:",paste(as.character(formulas[[i]])[-2],collapse=" "),sep=" "),
         xlab="coefficient",ylab="FE and 95% CI")
    segments(1:ncoef,mod_glm[[i]]$coefficients-1.96*se,1:ncoef,mod_glm[[i]]$coefficients+1.96*se)
    axis(1,at=1:ncoef,labels=names(coef(mod_glm[[i]])))
}
dev.off()

pdf(paste0(savedir,"graphs/fe_re_glms/all_cause_poisson_glm_residual_plots.pdf"),width=16,height=9)
for (i in 1:length(mod_glm)) {
    for (j in dims) {
        if (j %in% c("reg","res")) {
            plot(mod_glm[[i]]$residuals~chn_all[,factor(get(j))],
                 main=paste("Poisson GLM:",paste(as.character(formulas[[i]])[-2],collapse=" "),sep=" "),
                 xlab=j,ylab="Residuals"
            )
        } else {
            plot(mod_glm[[i]]$residuals~chn_all[,get(j)],
                 main=paste("Poisson GLM:",paste(as.character(formulas[[i]])[-2],collapse=" "),sep=" "),
                 xlab=j,ylab="Residuals")   
        }
    }
}
dev.off()

pdf(paste0(savedir,"graphs/fe_re_glms/all_cause_nb_glm_fe_plots.pdf"),width=16,height=9)
for (i in 1:length(mod_glm_nb)) {
    ncoef <- length(coef(mod_glm_nb[[i]]))
    se <- sqrt(diag(vcov(mod_glm_nb[[i]])))
    plot(coef(mod_glm_nb[[i]]),
         ylim=range(c(coef(mod_glm_nb[[i]])-1.96*se,coef(mod_glm_nb[[i]])+1.96*se)),
         xaxt="n",
         main=paste("NegBin GLM:",paste(as.character(formulas[[i]])[-2],collapse=" "),sep=" "),
         xlab="coefficient",ylab="FE and 95% CI")
    segments(1:ncoef,coef(mod_glm_nb[[i]])-1.96*se,1:ncoef,coef(mod_glm_nb[[i]])+1.96*se)
    axis(1,at=1:ncoef,labels=names(coef(mod_glm_nb[[i]])))
}
dev.off()

pdf(paste0(savedir,"graphs/fe_re_glms/all_cause_nb_glm_residual_plots.pdf"),width=16,height=9)
for (i in 1:length(mod_glm_nb)) {
    X <- model.matrix(formulas[[i]], data = chn_all)
    preds <- X %*% coef(mod_glm_nb[[i]])
    resids <- preds - chn_all$logmx
    plot(resids,
         main=paste("NegBin GLM:",paste(as.character(formulas[[i]])[-2],collapse=" "),sep=" "),
         ylab="Residuals",xlab="Observation")
    for (j in dims) {
        plot(resids~chn_all[,factor(get(j))],
             main=paste("NegBin GLM:",paste(as.character(formulas[[i]])[-2],collapse=" "),sep=" "),
             xlab=j,ylab="Residuals")
    }
}
dev.off()

## GLMMs ####

mod_glmm <- list()
mod_glmm_nb <- list()
formulas <- list()
i <- 0
for (r in 1:(length(loops)-1)) {
    for (rr in 1:nrow(loops[[r]])) {
        i <- i+1
        formulas[[i]] <- as.formula(paste("deaths ~ ",
                              paste(loops[[r]][rr,],collapse=" + "),
                              " + ", 
                              paste(paste("(1 | ",dims[which(!(dims %in% loops[[r]][rr,]))],")",sep=""),collapse=" + "),
                              sep=""))
        mod_glmm[[i]] <- glmer(formulas[[i]],offset=log(chn_all$exposure),data=chn_all,family=poisson(link=log))
        mod_glmm_nb[[i]] <- glmmTMB(formulas[[i]],offset=log(chn_all$exposure),data=chn_all,family=nbinom1(link="log"))
    }
}

## which has lowest AIC?
lapply(mod_glmm_nb,AIC)
best.glmm.nb <- mod_glmm_nb[[which(unlist(lapply(mod_glmm_nb,AIC))==min(unlist(lapply(mod_glmm_nb,AIC)),na.rm=T))]]

## compare glm and glmm AIC
AIC(best.glm.nb,best.glmm.nb)

## Random slope on year
chn_all$id <- paste(chn_all$agegp_name,chn_all$reg,chn_all$res,sep="_")
mod.glm.test1 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (year | id),offset=log(chn_all$exposure),data=chn_all,family=nbinom1(link="log"))
mod.glm.test2 <- glmmTMB(deaths ~ agegp_name + reg + res + (year | id),offset=log(chn_all$exposure),data=chn_all,family=nbinom1(link="log"))
mod.glm.test3 <- glmmTMB(deaths ~ agegp_name + reg + res + year + ar1(-1 + factor(year) | id),offset=log(chn_all$exposure),data=chn_all,family=nbinom1(link="log"))

AIC(mod.glm.test1,mod.glm.test2,best.glm.nb,best.glmm.nb,mod.glm.test3)

## additional models
mod3 <- glmmTMB(deaths ~ agegp_name + reg + res + (1 | year),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod5 <- glmmTMB(deaths ~ agegp_name + reg + res + year,offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod6 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (1 | year),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod7 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (1 | year) + (1 | agegp_name),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod8 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (1 | agegp_name),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod9 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (1 | agegp_name) + (1 | year) + (1 | reg) + (1 | res),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))

## plot residuals, FEs, and REs

pdf(paste0(savedir,"graphs/fe_re_glms/all_cause_nb_glmm_select_residual_comparisons.pdf"),width=16,height=9)

# plot predictions with REs
X <- model.matrix(~ agegp_name + reg + res,data=chn_all)
REs <- cbind(year=unique(chn_all$year),RE=ranef(mod3)$cond$year)
names(REs) <- c("year","RE")
results <- merge(chn_all,REs,by="year")
results$strata <- paste(results$res,results$reg,sep="_")
results <- results[order(results$strata,results$agegp,results$year),]
preds <- (X %*% fixef(mod3)$cond + results$RE)+log(results$exposure)
plot(preds-chn_all$logmx,ylim=c(-3,5))
abline(0,0,col="red")
plot((X %*% fixef(mod3)$cond)-chn_all$logmx,ylim=c(-3,5),
     main=paste("NegBin GLM:",paste(as.character(mod3$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
abline(0,0,col="red")
plot(preds~chn_all$logmx,main=paste("NegBin GLM:",paste(as.character(mod3$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))

X <- model.matrix(~ agegp_name + reg + res + year,data=chn_all)
preds3 <- X %*% fixef(mod5)$cond
plot(preds3-chn_all$logmx,ylim=c(-3,5),
     main=paste("NegBin GLM:",paste(as.character(mod5$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
abline(0,0,col="red")

plot(preds3~preds, ylab="No RE, FE on year, age, reg, res",xlab="RE on year, FE on age, reg, res",
     main="Prediction comparison")
abline(0,1,col="red")

# plot(preds3~preds2, ylab="No RE, FE on year, age, reg, res",xlab="RE on year and age, FE on year, age, reg, res",
#      main="Prediction comparison")
# abline(0,1,col="red")

# plot(preds2-preds)
# plot(preds3-preds)
# plot(preds3-preds2)

plot((preds3[chn_all$logmx!=-Inf]-chn_all$logmx[chn_all$logmx!=-Inf])~chn_all$agegp_name[chn_all$logmx!=-Inf],
     main="No RE, FE on year, age, reg, res")
plot((preds3[chn_all$logmx!=-Inf]-chn_all$logmx[chn_all$logmx!=-Inf])~jitter(chn_all$year[chn_all$logmx!=-Inf],1),
     main="No RE, FE on year, age, reg, res")
plot((preds3[chn_all$logmx!=-Inf]-chn_all$logmx[chn_all$logmx!=-Inf])~factor(chn_all$reg)[chn_all$logmx!=-Inf],
     main="No RE, FE on year, age, reg, res")
plot((preds3[chn_all$logmx!=-Inf]-chn_all$logmx[chn_all$logmx!=-Inf])~factor(chn_all$res)[chn_all$logmx!=-Inf],
     main="No RE, FE on year, age, reg, res")

plot((preds[chn_all$logmx!=-Inf]-chn_all$logmx[chn_all$logmx!=-Inf])~chn_all$agegp_name[chn_all$logmx!=-Inf],
     main="RE on year, FE on age, reg, res")
plot((preds[chn_all$logmx!=-Inf]-chn_all$logmx[chn_all$logmx!=-Inf])~jitter(chn_all$year[chn_all$logmx!=-Inf],1),
     main="RE on year, FE on age, reg, res")
plot((preds[chn_all$logmx!=-Inf]-chn_all$logmx[chn_all$logmx!=-Inf])~factor(chn_all$reg)[chn_all$logmx!=-Inf],
     main="RE on year, FE on age, reg, res")
plot((preds[chn_all$logmx!=-Inf]-chn_all$logmx[chn_all$logmx!=-Inf])~factor(chn_all$res)[chn_all$logmx!=-Inf],
     main="RE on year, FE on age, reg, res")

dev.off()

## plot GLMM results

pdf(paste0(savedir,"graphs/fe_re_glms/all_cause_nb_glmm_fe_plots.pdf"),width=16,height=9)
for (i in 1:length(mod_glmm_nb)) {
    ncoef <- length(coef(mod_glmm_nb[[i]]))
    se <- sqrt(diag(vcov(mod_glmm_nb[[i]])[[1]]))
    plot(fixef(mod_glmm_nb[[i]])$cond,
         ylim=range(c(fixef(mod_glmm_nb[[i]])$cond-1.96*se,fixef(mod_glmm_nb[[i]])$cond+1.96*se)),
         xaxt="n",main=paste("Poisson GLM:",paste(as.character(mod_glmm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "),
         xlab="coefficient",ylab="FE and 95% CI")
    segments(1:ncoef,fixef(mod_glmm_nb[[i]])$cond-1.96*se-1.96*se,1:ncoef,fixef(mod_glmm_nb[[i]])$cond-1.96*se+1.96*se)
    axis(1,at=1:ncoef,labels=names(coef(mod_glmm_nb[[i]])))
}
dev.off()

# Cause-specific mortality analyses ####

## GLMs ####
dims2 <- c("year","agegp_name","reg","res","cause")

loops2 <- list()
for (rr in 1:length(dims2)) {
    loops2[[rr]] <- combinations(n=length(dims2),r=rr,v=dims2)
}

mod_cause_glm <- list()
formulals <- list()
i <- 0
for (r in 1:length(loops2)) {
    for (rr in 1:nrow(loops2[[r]])) {
        i <- i+1
        formulas[[i]] <- as.formula(paste("deaths ~ offset(logpy) + ",paste(loops2[[r]][rr,],collapse=" + "),sep=""))
        mod_cause_glm[[i]] <- glm(formulas[[i]],data=chn_all, family = poisson)
    }
}

## AIC analysis
do.call(anova, mod_cause_glm)
best.cause.glm <- mod_cause_glm[[which(unlist(lapply(mod_cause_glm,AIC))==min(unlist(lapply(mod_cause_glm,AIC))))]]

## plot FEs and residuals
pdf(paste0(savedir,"graphs/fe_re_glms/cause_specific_glm_fe_plots.pdf"),width=16,height=9)
par(oma=c(5,2,2,2),xpd=NA)
for (i in 1:length(mod_cause_glm)) {
    ncoef <- length(coef(mod_cause_glm[[i]]))
    se <- sqrt(diag(vcov(mod_cause_glm[[i]])))
    plot(coef(mod_cause_glm[[i]]),
         ylim=range(c(coef(mod_cause_glm[[i]])-1.96*se,coef(mod_cause_glm[[i]])+1.96*se)),
         xaxt="n",
         main=paste("Pois GLM:",paste(as.character(formulas[[i]])[-2],collapse=" "),sep=" "),
         xlab="",ylab="FE and 95% CI")
    segments(1:ncoef,coef(mod_cause_glm[[i]])-1.96*se,1:ncoef,coef(mod_cause_glm[[i]])+1.96*se)
    axis(1,at=1:ncoef,labels=FALSE)
    text(x=1:ncoef,y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),labels=names(coef(mod_cause_glm[[i]])),
         srt=45, adj=1, xpd=NA)
}
dev.off()

pdf(paste0(savedir,"graphs/fe_re_glms/cause_specific_glm_residual_plots_mx.pdf"),width=16,height=9)
for (i in 1:length(mod_cause_glm)) {
    X <- model.matrix(formulas[[i]],data=chn_all)
    preds <- exp(X %*% coef(mod_cause_glm[[i]]))
    resids <- preds - chn_all$mx
    plot(resids,ylab="Residuals (mx scale)",xlab="Observation",
         main=paste("Pois GLM:",paste(as.character(formulas[[i]])[-2],collapse=" "),sep=" "))
    for (j in dims2) {
        plot(resids~chn_all[,factor(get(j))],xlab=j,ylab="Residuals",
             main=paste("Pois GLM:",paste(as.character(formulas[[i]])[-2],collapse=" "),sep=" "))
    }
}
dev.off()

pdf(paste0(savedir,"graphs/fe_re_glms/cause_specific_glm_residual_plots_logmx_noInf.pdf"),width=16,height=9)
for (i in 1:length(mod_cause_glm)) {
    X <- model.matrix(formulas[[i]],data=chn_all)
    preds <- X %*% coef(mod_cause_glm[[i]])
    resids <- preds[!is.infinite(chn_all$logmx)] - chn_all$logmx[!is.infinite(chn_all$logmx)]
    plot(resids,ylab="Residuals (log scale, no -Inf)",xlab="Observation",
         main=paste("Pois GLM:",paste(as.character(formulas[[i]])[-2],collapse=" "),sep=" "))
    for (j in dims2) {
        plot(resids~chn_all[!is.infinite(chn_all$logmx),factor(get(j))],xlab=j,ylab="Residuals",
             main=paste("Pois GLM:",paste(as.character(formulas[[i]])[-2],collapse=" "),sep=" "))
    }
}
dev.off()

## cause correlations in residuals
big_i <- length(mod_cause_glm)
X <- model.matrix(formulas[[big_i]],data=chn_all)
preds <- X %*% coef(mod_cause_glm[[big_i]])
chn_all$preds <- preds
chn_all$exp_resids <- exp(chn_all$preds) - chn_all$mx
chn_all$log_resids <- chn_all$preds - chn_all$logmx

cause_resids <- as.data.frame(reshape(chn_all[chn_all$logmx!=-Inf,c("log_resids","exp_resids","cause","agegp_name","year","reg","res")],
                                      direction="wide",v.names=c("log_resids","exp_resids"),timevar=c("cause"),
                                      idvar=c("agegp_name","year","reg","res")))
cor(cause_resids[,grep("exp_resids",names(cause_resids),value=T)],use="pairwise.complete.obs")
cor(cause_resids[,grep("log_resids",names(cause_resids),value=T)],use="pairwise.complete.obs")

deathvars <- unique(chn_all$cause)
pdf(paste0(savedir,"graphs/fe_re_glms/cause_resid_pairs_plot.pdf"),width=16,height = 9)
for (i in (1:length(causes))) {
    for (j in 1:i) {
        if (i == j) next
        c1 <- deathvars[i]
        c2 <- deathvars[j]
        plot(cause_resids[,paste0("log_resids.",c1)] ~ cause_resids[,paste0("log_resids.",c2)],
             xlab=c2,ylab=c1)
    }
}
pairs(cause_resids[,grep("log_resids",names(cause_resids),value=T)])
dev.off()

## GLMMs ####

mod_cause_glmm_nb <- list()

i <- 0
for (r in 1:(length(loops2)-1)) {
    for (rr in 1:nrow(loops2[[r]])) {
        i <- i+1
        f <- as.formula(paste("deaths ~ offset(logpy) +",
                              paste(loops2[[r]][rr,],collapse=" + "),
                              " + ", 
                              paste(paste("(1 | ",dims2[which(!(dims2 %in% loops2[[r]][rr,]))],")",sep=""),collapse=" + "),
                              sep=""))
        mod_cause_glmm_nb[[i]] <- glmmTMB(f,offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
    }
}

## which has lowest AIC?
lapply(mod_cause_glmm_nb,AIC)
best.cause.glmm.nb <- mod_cause_glmm_nb[[which(unlist(lapply(mod_cause_glmm_nb,AIC))==min(unlist(lapply(mod_cause_glmm_nb,AIC))))]]

## compare glm and glmm AIC
AIC(best.cause.glm,best.cause.glmm.nb)

## extra models
# chn_all$id <- paste(chn_all$agegp_name,chn_all$reg,chn_all$res,chn_all$cause,sep="_")
# mod1 <- glmmTMB(deaths ~ agegp_name + cause + (year | id),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
# mod2 <- glmmTMB(deaths ~ agegp_name + cause + reg + res + (year | id),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
# mod3 <- glmmTMB(deaths ~ agegp_name + cause + reg + res + year + (year | id),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
# mod4 <- glmmTMB(deaths ~ agegp_name + cause + reg + res + year + (1 + year | id),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
# 
# AIC(mod1,mod2,mod3,mod4)
# 
# chn_all$id2 <- paste(chn_all$agegp_name,chn_all$reg,chn_all$res,chn_all$year,sep="_")
# mod4 <- glmmTMB(deaths ~ agegp_name + cause + reg + res + year + (cause | id2),offset=log(chn_all$exposure),data=chn_all,family=nbinom1(link="log"))
# mod5 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (cause | id2),offset=log(chn_all$exposure),data=chn_all,family=nbinom1(link="log"))
# mod6 <- glmmTMB(deaths ~ agegp_name + cause + reg + res + year + (cause | id2) + (year | id),offset=log(chn_all$exposure),data=chn_all,family=nbinom1(link="log"))
# mod7 <- glmmTMB(deaths ~ agegp_name + reg + res + (cause | id2) + (year | id),offset=log(chn_all$exposure),data=chn_all,family=nbinom1(link="log"))
# mod8 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (cause | id2) + (year | id),offset=log(chn_all$exposure),data=chn_all,family=nbinom1(link="log"))
# mod9 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (cause | id2),offset=log(chn_all$exposure),data=chn_all,family=nbinom1(link="log"))

mod0 <- glmmTMB(deaths ~ agegp_name + reg + res + year + cause,
                offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod1 <- glmmTMB(deaths ~ (1 | year) + (1 | agegp) + (1 | reg) + (1 | res) + (1 | cause),
                offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
chn_all$id <- paste(chn_all$agegp_name,chn_all$reg,chn_all$res,chn_all$year,sep="_")
mod2 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (-1 + cause | id),
                offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
chn_all$id2 <- paste(chn_all$agegp_name,chn_all$reg,chn_all$res,sep="_")
mod2a <- glmmTMB(deaths ~ agegp_name + reg + res + year + (-1 + cause | id2),
                 offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
anova(mod2,mod2a)
ranef(mod_cause_glmm_nb[[29]])
ranef(mod2)

anova(mod2,mod_cause_glmm_nb[[29]])
plot(residuals(mod2) ~ factor(chn_all$year))
plot(residuals(mod2) ~ factor(chn_all$agegp_name))
plot(residuals(mod2) ~ factor(chn_all$reg))
plot(residuals(mod2) ~ factor(chn_all$res))
plot(residuals(mod2) ~ factor(chn_all$cause))

chn_all$id2 <- paste(chn_all$agegp_name,chn_all$reg,chn_all$res,sep="_")
mod3 <- glmmTMB(deaths ~ agegp_name + reg + res + (1 + year + cause | id2),
                offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
anova(mod2,mod3)

## cause correlations in residuals
results <- chn_all
results$preds <- predict(mod3)
results$resids <- residuals(mod3)

cause_resids2 <- as.data.frame(reshape(results[,c("resids","preds","cause","agegp_name","year","reg","res")],
                                       direction="wide",v.names=c("resids","preds"),timevar=c("cause"),
                                       idvar=c("agegp_name","year","reg","res")))
cor(cause_resids2[,grep("resids",names(cause_resids2),value=T)],use="pairwise.complete.obs")
cor(cause_resids2[,grep("preds",names(cause_resids2),value=T)],use="pairwise.complete.obs")

cause_resids3 <- as.data.frame(reshape(results[,c("resids","preds","cause","agegp_name","year","reg","res")],
                                       direction="wide",v.names=c("resids","preds"),timevar=c("agegp_name"),
                                       idvar=c("cause","year","reg","res")))
cor(cause_resids3[,grep("resids",names(cause_resids3),value=T)],use="pairwise.complete.obs")
cor(cause_resids3[,grep("preds",names(cause_resids3),value=T)],use="pairwise.complete.obs")

pairs(cause_resids3[,grep("resids",names(cause_resids3),value=T)],
      main=paste("NegBin GLM:",paste(as.character(mod3$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))

pdf(paste0(savedir,"graphs/fe_re_glms/cause_resid_pairs_plot_glmm_yearRE_corcauseRE.pdf"),width=16,height = 9)
pairs(cause_resids2[,grep("resids",names(cause_resids2),value=T)],
      main=paste("NegBin GLM:",paste(as.character(mod3$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
dev.off()

# Testing different dimensions for time RE ####

## raw data
chn_agg_age_year <- aggregate(cbind(chn_all$deaths,chn_all$exposure),
                              by=list(chn_all$agegp,chn_all$year),sum)
chn_agg_age_year$mx <- chn_agg_age_year[,3]/chn_agg_age_year[,4]
chn_agg_age_year$col <- as.numeric(as.factor(chn_agg_age_year$Group.1))

chn_agg_cause_year <- aggregate(cbind(chn_all$deaths,chn_all$exposure),
                                by=list(chn_all$cause,chn_all$year),sum)
chn_agg_cause_year$mx <- chn_agg_cause_year[,3]/chn_agg_cause_year[,4]
chn_agg_cause_year$col <- as.numeric(as.factor(chn_agg_cause_year$Group.1))

chn_all$strata <- paste(chn_all$reg, chn_all$res, sep = "_")
chn_agg_strata_year <- aggregate(cbind(chn_all$deaths,chn_all$exposure),
                                 by=list(chn_all$strata,chn_all$year),sum)
chn_agg_strata_year$mx <- chn_agg_strata_year[,3]/chn_agg_strata_year[,4]
chn_agg_strata_year$col <- as.numeric(as.factor(chn_agg_strata_year$Group.1))

plot(mx~Group.2,data=chn_agg_age_year,col=chn_agg_age_year$col)
plot(mx~Group.2,data=chn_agg_cause_year,col=chn_agg_cause_year$col)
plot(mx~Group.2,data=chn_agg_strata_year,col=chn_agg_strata_year$col)

## modeling
chn_all$ageXcause <- paste(chn_all$agegp,chn_all$cause,sep="_")
chn_all$ageXcause <- ifelse(chn_all$agegp==1 | chn_all$cause=="nCH10" | (chn_all$agegp %in% c(4,5,6) & chn_all$cause=="nCH17"), "ref",chn_all$ageXcause)
chn_all$ageXcause <- factor(chn_all$ageXcause)
chn_all$ageXcause <- relevel(chn_all$ageXcause,ref="ref")
chn_all$id1 <- paste(chn_all$agegp_name,chn_all$reg,chn_all$res,sep="_")
chn_all$id2 <- paste(chn_all$agegp_name,chn_all$cause,sep="_")
chn_all$id3 <- paste(chn_all$reg,chn_all$res,chn_all$cause,sep="_")
chn_all$id4 <- paste(chn_all$agegp_name,chn_all$cause,chn_all$reg,chn_all$res,sep="_")
chn_all$id5 <- paste(chn_all$agegp_name,chn_all$year,chn_all$reg,chn_all$res,sep="_")

mod1 <- glmmTMB(deaths ~ year + agegp_name + cause + reg + res + (1 + year|agegp_name),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod2 <- glmmTMB(deaths ~ year + agegp_name + cause + reg + res + (1 + year|cause),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod3 <- glmmTMB(deaths ~ year + agegp_name + cause + reg + res + (1 + year|strata),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod4 <- glmmTMB(deaths ~ year + agegp_name + cause + reg + res + (1 + year|id1),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod5 <- glmmTMB(deaths ~ year + agegp_name + cause + reg + res + (1 + year|id2),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod6 <- glmmTMB(deaths ~ year + agegp_name + cause + reg + res + (1 + year|id3),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))
mod7 <- glmmTMB(deaths ~ year + agegp_name + cause + reg + res + (1 + year|id4),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))

mod8 <- glmmTMB(deaths ~ year + agegp_name + cause + ageXcause + reg + res + (1 + year|id1),offset=log(chn_all$exposure),data=chn_all,family=nbinom2(link="log"))

anova(mod1,mod2,mod3,mod4,mod5,mod6,mod7)
anova(mod5,mod7)
anova(mod7,mod8)

res <- list()
length(res) <- 7
for(i in 1:7) {
    chn_all[,paste0("res_mod",i)] <- residuals(get(paste0("mod",i)))
}

for (i in 1:7) {
    print(xyplot(chn_all[,get(paste0("res_mod",i))]~chn_all$year | chn_all$agegp_name,groups=chn_all$cause,main=paste0("mod",i,": colors = cause")))
    print(xyplot(chn_all[,get(paste0("res_mod",i))]~chn_all$year | chn_all$agegp_name,groups=chn_all$strata,main=paste0("mod",i,": colors = strata")))
    print(xyplot(chn_all[,get(paste0("res_mod",i))]~chn_all$year | chn_all$cause,groups=chn_all$agegp_name,main=paste0("mod",i,": colors = age")))
    print(xyplot(chn_all[,get(paste0("res_mod",i))]~chn_all$year | chn_all$cause,groups=chn_all$strata,main=paste0("mod",i,": colors = strata")))
    print(xyplot(chn_all[,get(paste0("res_mod",i))]~chn_all$year | factor(chn_all$strata),groups=chn_all$agegp_name,main=paste0("mod",i,": colors = age")))
    print(xyplot(chn_all[,get(paste0("res_mod",i))]~chn_all$year | factor(chn_all$strata),groups=chn_all$cause,main=paste0("mod",i,": colors = cause")))
}

# Plot model predictions and compare to raw data ####
causes <- unique(chn_all$cause)
stratas <- unique(chn_all$strata)
ages <- unique(chn_all$agegp)
chn_all$preds <- predict(mod7,type="link") - log(chn_all$exposure)

pdf(paste0(savedir,"graphs/fe_re_glms/glmm_model_preds_vs_data_over_time_noagecauseinteraction.pdf"),width=16,height=9)
for (i in 1:length(unique(chn_all$strata))) {
    for (j in 1:length(unique(chn_all$cause))) {
        par(mfrow=c(2,3))
        for (k in 1:length(unique(chn_all$agegp))) {
            # skip impossible age-cause combos
            if (ages[k] %in% c(4, 5, 6) & causes[j] %in% c("nCH10", "nCH11")) next
            
            # subset data and plot
            tmp <- chn_all[chn_all$strata==stratas[i] & chn_all$cause == causes[j] & chn_all$agegp == ages[k], ]
            plot(logmx~year,data=tmp,col=alpha("red",0.5),pch=19,ylim=range(c(tmp$preds,tmp$logmx[!is.infinite(tmp$logmx)])),
                 main=paste0("cause: ",causes[j],"; strata: ",stratas[i],"; age group: ",ages[k]))
            lines(preds~year,data=tmp,col=alpha("blue",0.5),pch=19)
        }
    }
}
dev.off() 
