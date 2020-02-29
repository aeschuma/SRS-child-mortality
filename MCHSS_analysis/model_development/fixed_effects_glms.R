## Austin Schumacher
## 2018 Nov 23
## 1. Fit GLMs on China data to assess the drivers of variability and structure/correlations in residuals
## 2. Fit GLMMs on China data to assess model fit, drivers of variability, and correlation in residuals

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

library(data.table);library(gtools);library(RColorBrewer);library(MASS);library(lme4);
library(glmmTMB);library(stringr);library(lattice);library(scales);

## define directories

# working directory for code
wd <- paste(root,"Desktop/dissertation/model_development/glm_glmm",sep="")

# directory to load MCHSS data
datadir <- paste(root,"Dropbox/dissertation_2/cause_specific_child_mort/china_explore",sep="")

## set directory
setwd(datadir)

## read in data
load("../china_explore/chn_srs_formatted.Rdata")

#####  format data
deathvars <- grep("d_n",names(chn_srs),value=T)

## list of causes
causes <- substr(deathvars,3,nchar(deathvars))

## names for causes
cause_names <- data.frame(cause=causes,
                          cause_name=c("prematurity","birth asphyxia/trauma",
                                       "congenital anomalies","other non-communicable",
                                       "injuries","diarrhea",
                                       "acute resp. infections","other grp 1"))

## data frame for old/young cause reference
cause_ref <- data.frame(cause=causes,
                        cause_name=c("prematurity","birth asphyxia/trauma",
                                     "congenital anomalies","other non-communicable",
                                     "injuries","diarrhea",
                                     "acute resp. infections","other grp 1"),
                        cause_ind=1:8,
                        cause_young_ind=1:8,
                        cause_old_ind=c(0,0,1:6))
causes_young <- causes
causes_old <- causes_young[!(causes_young %in% c("nCH10","nCH11"))]
chn_srs$strata <- paste(chn_srs$res,chn_srs$reg,sep="_")
ages <- unique(chn_srs$agegp)
young_ages <- 1:3
old_ages <- 4:6
years <- unique(chn_srs$year)
stratas <- unique(chn_srs$strata)

## format data long
chn <- reshape(chn_srs,direction="long",varying=deathvars,v.names="deaths",timevar="cause",
               times=substr(deathvars,3,length(deathvars)),
               idvar=c("year","agegp","reg","res"))
row.names(chn) <- NULL

## extra formatting
chn$deaths <- ceiling(chn$deaths)
chn$mx <- chn$deaths/chn$exposure
chn$logmx <- log(chn$mx)
chn <- chn[order(chn$strata,chn$agegp,chn$year,chn$cause),]
chn_all <- as.data.table(chn[!((chn$agegp %in% c(4,5,6)) & (chn$cause %in% c("nCH10","nCH11"))),])

chn <- chn_all[,list(exposure=mean(exposure),
                     deaths=sum(deaths),
                     mx=sum(deaths)/mean(exposure),
                     logmx=log(sum(deaths)/mean(exposure))),
               by=list(agegp,year,reg,res)]
chn$agegp_name <- factor(chn$agegp)
chn$logpy <- log(chn$exposure)

chn_all <- as.data.frame(chn_all)
chn_all$logpy <- log(chn_all$exposure)
chn_all$young_age <- ifelse(chn_all$agegp %in% 1:2,1,0)
chn_all <- merge(chn_all,cause_ref)
chn_all$prematurity <- ifelse(chn_all$cause_name=="prematurity",1,0)
chn_all$ageXcause <- paste(chn_all$agegp,chn_all$cause,sep="_")
chn_all$ageXcause <- ifelse(chn_all$agegp==1 | chn_all$cause=="nCH10" | (chn_all$agegp %in% c(4,5,6) & chn_all$cause=="nCH17"), "ref",chn_all$ageXcause)
chn_all$ageXcause <- factor(chn_all$ageXcause)
chn_all$ageXcause <- relevel(chn_all$ageXcause,ref="ref")
chn_all <- chn_all[order(chn_all$reg,chn_all$res,chn_all$agegp,chn_all$year,chn_all$cause_ind),]

## GLMs on all cause mortality
dims <- c("year","agegp_name","reg","res")

loops <- list()
for (rr in 1:length(dims)) {
    loops[[rr]] <- combinations(n=length(dims),r=rr,v=dims)
}
mod_glm <- list()
mod_glm_nb <- list()

i <- 0
for (r in 1:length(loops)) {
    for (rr in 1:nrow(loops[[r]])) {
        i <- i+1
        f <- as.formula(paste("deaths ~ offset(logpy) + ",paste(loops[[r]][rr,],collapse=" + "),sep=""))
        mod_glm[[i]] <- glm(f,data=chn,family=quasipoisson(link="log"))
        mod_glm_nb[[i]] <- glm.nb(f,data=chn)
    }
}

# look at residuals
plot(mod_glm[[1]]$residuals~chn$year)
plot(mod_glm[[1]]$residuals~chn$agegp_name)
plot(mod_glm[[4]]$residuals~chn$year)
plot(mod_glm[[4]]$residuals~chn$agegp_name)
plot(mod_glm[[5]]$residuals~chn$year)
plot(mod_glm[[5]]$residuals~chn$agegp_name)

cor(mod_glm[[11]]$residuals,chn$year,method="spearman")
cor(mod_glm[[14]]$residuals,as.numeric(chn$agegp_name),method="spearman")
cor(mod_glm[[12]]$residuals,as.numeric(factor(chn$res)),method="spearman")
cor(mod_glm[[13]]$residuals,as.numeric(factor(chn$reg)),method="spearman")

# look at residuals
# plot(predict(mod_glm_nb[[1]])-log(chn$deaths)~chn$year)
# plot(predict(mod_glm_nb[[1]])-log(chn$deaths)~chn$agegp_name)
# plot(predict(mod_glm_nb[[4]])-log(chn$deaths)~chn$year)
# plot(predict(mod_glm_nb[[4]])-log(chn$deaths)~chn$agegp_name)
# plot(predict(mod_glm_nb[[5]])-log(chn$deaths)~chn$year)
# plot(predict(mod_glm_nb[[5]])-log(chn$deaths)~chn$agegp_name)
# 
# cor(predict(mod_glm_nb[[11]])-log(chn$deaths),chn$year,method="spearman")
# cor(predict(mod_glm_nb[[14]])-log(chn$deaths),as.numeric(chn$agegp_name),method="spearman")
# cor(predict(mod_glm_nb[[12]])-log(chn$deaths),as.numeric(factor(chn$res)),method="spearman")
# cor(predict(mod_glm_nb[[13]])-log(chn$deaths),as.numeric(factor(chn$reg)),method="spearman")
# 
# cor(cbind((predict(mod_glm_nb[[13]])-log(chn$deaths))[chn$reg=="east"],
#           (predict(mod_glm_nb[[13]])-log(chn$deaths))[chn$reg=="mid"],
#           (predict(mod_glm_nb[[13]])-log(chn$deaths))[chn$reg=="west"]),
#     method="spearman")

## AIC analysis
lapply(mod_glm_nb,AIC)
best.glm.nb <- mod_glm_nb[[which(unlist(lapply(mod_glm_nb,AIC))==min(unlist(lapply(mod_glm_nb,AIC))))]]

## plot FEs and residuals
setwd(wd)
pdf("graphs/all_cause_poisson_glm_fe_plots.pdf",width=16,height=9)
for (i in 1:length(mod_glm)) {
    ncoef <- length(mod_glm[[i]]$coefficients)
    se <- summary(mod_glm[[i]])$coef[,2]
    plot(mod_glm[[i]]$coefficients,ylim=range(c(mod_glm[[i]]$coefficients-1.96*se,mod_glm[[i]]$coefficients+1.96*se)),
         xaxt="n",main=paste("Poisson GLM:",paste(as.character(mod_glm[[i]]$formula)[-2],collapse=" "),sep=" "),
         xlab="coefficient",ylab="FE and 95% CI")
    segments(1:ncoef,mod_glm[[i]]$coefficients-1.96*se,1:ncoef,mod_glm[[i]]$coefficients+1.96*se)
    axis(1,at=1:ncoef,labels=names(coef(mod_glm[[i]])))
}
dev.off()

pdf("graphs/all_cause_poisson_glm_residual_plots.pdf",width=16,height=9)
for (i in 1:length(mod_glm)) {
    for (j in dims) {
        if (j %in% c("reg","res")) {
            plot(mod_glm[[i]]$residuals~chn[,factor(get(j))],xlab=j,ylab="Residuals",
                 main=paste("Poisson GLM:",paste(as.character(mod_glm[[i]]$formula)[-2],collapse=" "),sep=" "))
        } else {
            plot(mod_glm[[i]]$residuals~chn[,get(j)],xlab=j,ylab="Residuals",
                 main=paste("Poisson GLM:",paste(as.character(mod_glm[[i]]$formula)[-2],collapse=" "),sep=" "))   
        }
    }
}
dev.off()

pdf("graphs/all_cause_nb_glm_fe_plots.pdf",width=16,height=9)
for (i in 1:length(mod_glm_nb)) {
    ncoef <- length(fixef(mod_glm_nb[[i]])$cond)
    se <- sqrt(diag(vcov(mod_glm_nb[[i]])$cond))
    plot(fixef(mod_glm_nb[[i]])$cond,ylim=range(c(fixef(mod_glm_nb[[i]])$cond-1.96*se,fixef(mod_glm_nb[[i]])$cond+1.96*se)),
         xaxt="n",main=paste("Poisson GLM:",paste(as.character(mod_glm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "),
         xlab="coefficient",ylab="FE and 95% CI")
    segments(1:ncoef,fixef(mod_glm_nb[[i]])$cond-1.96*se,1:ncoef,fixef(mod_glm_nb[[i]])$cond+1.96*se)
    axis(1,at=1:ncoef,labels=names(coef(mod_glm_nb[[i]])))
}
dev.off()

pdf("graphs/all_cause_nb_glm_residual_plots.pdf",width=16,height=9)
for (i in 1:length(mod_glm_nb)) {
    X <- model.matrix(mod_glm_nb[[i]]$modelInfo$allForm$formula,data=chn)
    preds <- X %*% fixef(mod_glm_nb[[i]])$cond
    resids <- preds - chn$logmx
    plot(resids,ylab="Residuals",xlab="Observation",
         main=paste("NegBin GLM:",paste(as.character(mod_glm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
    for (j in dims) {
        plot(resids~chn[,factor(get(j))],xlab=j,ylab="Residuals",
             main=paste("NegBin GLM:",paste(as.character(mod_glm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
    }
}
dev.off()

## GLMMs on all-cause mortality

mod_glmm <- list()
mod_glmm_nb <- list()

i <- 0
for (r in 1:(length(loops)-1)) {
    for (rr in 1:nrow(loops[[r]])) {
        i <- i+1
        f <- as.formula(paste("deaths ~ ",
                              paste(loops[[r]][rr,],collapse=" + "),
                              " + ", 
                              paste(paste("(1 | ",dims[which(!(dims %in% loops[[r]][rr,]))],")",sep=""),collapse=" + "),
                              sep=""))
        mod_glmm[[i]] <- glmer(f,offset=log(chn$exposure),data=chn,family=poisson(link=log))
        mod_glmm_nb[[i]] <- glmmTMB(f,offset=log(chn$exposure),data=chn,family=nbinom1(link="log"))
    }
}

## which has lowest AIC?
lapply(mod_glmm_nb,AIC)
best.glmm.nb <- mod_glmm_nb[[which(unlist(lapply(mod_glmm_nb,AIC))==min(unlist(lapply(mod_glmm_nb,AIC)),na.rm=T))]]

## compare glm and glmm AIC
AIC(best.glm.nb,best.glmm.nb)

## Random slope on year
chn$id <- paste(chn$agegp_name,chn$reg,chn$res,sep="_")
mod.glm.test1 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (year | id),offset=log(chn$exposure),data=chn,family=nbinom1(link="log"))
mod.glm.test2 <- glmmTMB(deaths ~ agegp_name + reg + res + (year | id),offset=log(chn$exposure),data=chn,family=nbinom1(link="log"))
mod.glm.test3 <- glmmTMB(deaths ~ agegp_name + reg + res + year + ar1(year | id),offset=log(chn$exposure),data=chn,family=nbinom1(link="log"))

AIC(mod.glm.test1,mod.glm.test2,best.glm.nb,best.glmm.nb,mod.glm.test3)

## additional models
mod3 <- glmmTMB(deaths ~ agegp_name + reg + res + (1 | year),offset=log(chn$exposure),data=chn,family=nbinom2(link="log"))
mod5 <- glmmTMB(deaths ~ agegp_name + reg + res + year,offset=log(chn$exposure),data=chn,family=nbinom2(link="log"))
mod6 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (1 | year),offset=log(chn$exposure),data=chn,family=nbinom2(link="log"))
mod7 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (1 | year) + (1 | agegp_name),offset=log(chn$exposure),data=chn,family=nbinom2(link="log"))
mod8 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (1 | agegp_name),offset=log(chn$exposure),data=chn,family=nbinom2(link="log"))
mod9 <- glmmTMB(deaths ~ agegp_name + reg + res + year + (1 | agegp_name) + (1 | year) + (1 | reg) + (1 | res),offset=log(chn$exposure),data=chn,family=nbinom2(link="log"))

## plot residuals, FEs, and REs

pdf("graphs/all_cause_nb_glmm_select_residual_comparisons.pdf",width=16,height=9)

# plot predictions with REs
X <- model.matrix(~ agegp_name + reg + res,data=chn)
REs <- cbind(year=years,RE=ranef(mod3)$cond$year)
names(REs) <- c("year","RE")
results <- merge(chn,REs,by="year")
results$strata <- paste(results$res,results$reg,sep="_")
results <- results[order(results$strata,results$agegp,results$year),]
preds <- (X %*% fixef(mod3)$cond + results$RE)+log(results$exposure)
plot(preds-chn$logmx,ylim=c(-3,5))
abline(0,0,col="red")
plot((X %*% fixef(mod3)$cond)-chn$logmx,ylim=c(-3,5),
     main=paste("NegBin GLM:",paste(as.character(mod3$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
abline(0,0,col="red")
plot(preds~chn$logmx,main=paste("NegBin GLM:",paste(as.character(mod3$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))

X <- model.matrix(~ agegp_name + reg + res + year,data=chn)
preds3 <- X %*% fixef(mod5)$cond
plot(preds3-chn$logmx,ylim=c(-3,5),
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

plot((preds3[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~chn$agegp_name[chn$logmx!=-Inf],
     main="No RE, FE on year, age, reg, res")
plot((preds3[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~jitter(chn$year[chn$logmx!=-Inf],1),
     main="No RE, FE on year, age, reg, res")
plot((preds3[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~factor(chn$reg)[chn$logmx!=-Inf],
     main="No RE, FE on year, age, reg, res")
plot((preds3[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~factor(chn$res)[chn$logmx!=-Inf],
     main="No RE, FE on year, age, reg, res")

plot((preds[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~chn$agegp_name[chn$logmx!=-Inf],
     main="RE on year, FE on age, reg, res")
plot((preds[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~jitter(chn$year[chn$logmx!=-Inf],1),
     main="RE on year, FE on age, reg, res")
plot((preds[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~factor(chn$reg)[chn$logmx!=-Inf],
     main="RE on year, FE on age, reg, res")
plot((preds[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~factor(chn$res)[chn$logmx!=-Inf],
     main="RE on year, FE on age, reg, res")

# plot((preds3[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~chn$agegp_name[chn$logmx!=-Inf])
# plot((preds3[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~jitter(chn$year[chn$logmx!=-Inf],1))
# plot((preds3[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~factor(chn$reg)[chn$logmx!=-Inf])
# plot((preds3[chn$logmx!=-Inf]-chn$logmx[chn$logmx!=-Inf])~factor(chn$res)[chn$logmx!=-Inf])
dev.off()

## plot GLMM results

# pdf("graphs/all_cause_nb_glmm_fe_plots.pdf",width=16,height=9)
# for (i in 1:length(mod_glmm_nb)) {
#     ncoef <- length(fixef(mod_glmm_nb[[i]])$cond)
#     se <- sqrt(diag(vcov(mod_glmm_nb[[i]])$cond))
#     plot(fixef(mod_glmm_nb[[i]])$cond,ylim=range(c(fixef(mod_glmm_nb[[i]])$cond-1.96*se,fixef(mod_glmm_nb[[i]])$cond+1.96*se)),
#          xaxt="n",main=paste("Poisson GLM:",paste(as.character(mod_glmm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "),
#          xlab="coefficient",ylab="FE and 95% CI")
#     segments(1:ncoef,fixef(mod_glmm_nb[[i]])$cond-1.96*se-1.96*se,1:ncoef,fixef(mod_glmm_nb[[i]])$cond-1.96*se+1.96*se)
#     axis(1,at=1:ncoef,labels=names(coef(mod_glmm_nb[[i]])))
# }
# dev.off()
# 
# pdf("graphs/all_cause_nb_glmm_residual_plots.pdf",width=16,height=9)
# for (i in 1:length(mod_glmm_nb)) {
#     i <- 1
#     f <- substr(as.character(mod_glmm_nb[[i]]$modelInfo$allForm$formula)[3],
#                 grep("\\(",as.character(mod_glmm_nb[[i]]$modelInfo$allForm$formula)[3]),
#                 )
#     X <- model.matrix(mod_glmm_nb[[i]]$modelInfo$allForm$formula,data=chn)
#     REs <- ranef(mod_glmm_nb[[i]])$cond$year
#     preds <- X %*% fixef(mod_glmm_nb[[i]])$cond
#     resids <- preds - chn$logmx
#     plot(resids,ylab="Residuals",xlab="Observation",
#          main=paste("NegBin GLM:",paste(as.character(mod_glmm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
#     for (j in dims) {
#         j <- "year"
#         plot(resids~chn[,factor(get(j))],xlab=j,ylab="Residuals",
#              main=paste("NegBin GLM:",paste(as.character(mod_glm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
#     }
# }
# dev.off()

####################
#### Cause specific mortality
####################

dims2 <- c("year","agegp_name","reg","res","cause")

loops2 <- list()
for (rr in 1:length(dims2)) {
    loops2[[rr]] <- combinations(n=length(dims2),r=rr,v=dims2)
}

mod_cause_glm <- list()

i <- 0
for (r in 1:length(loops2)) {
    for (rr in 1:nrow(loops2[[r]])) {
        i <- i+1
        f <- as.formula(paste("deaths ~ offset(logpy) + ",paste(loops2[[r]][rr,],collapse=" + "),sep=""))
        mod_cause_glm[[i]] <- glm(f,data=chn_all,family=quasipoisson(link="log"))
    }
}

## AIC analysis
do.call(anova,mod_cause_glm)
best.cause.glm <- mod_cause_glm_nb[[which(unlist(lapply(mod_cause_glm_nb,AIC))==min(unlist(lapply(mod_cause_glm_nb,AIC))))]]

## plot FEs and residuals
setwd(wd)
pdf("graphs/cause_specific_nb_glm_fe_plots.pdf",width=16,height=9)
par(oma=c(5,2,2,2),xpd=NA)
for (i in 1:length(mod_cause_glm_nb)) {
    ncoef <- length(fixef(mod_cause_glm_nb[[i]])$cond)
    se <- sqrt(diag(vcov(mod_cause_glm_nb[[i]])$cond))
    plot(fixef(mod_cause_glm_nb[[i]])$cond,ylim=range(c(fixef(mod_cause_glm_nb[[i]])$cond-1.96*se,fixef(mod_cause_glm_nb[[i]])$cond+1.96*se)),
         xaxt="n",main=paste("NegBin GLM:",paste(as.character(mod_cause_glm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "),
         xlab="",ylab="FE and 95% CI")
    segments(1:ncoef,fixef(mod_cause_glm_nb[[i]])$cond-1.96*se,1:ncoef,fixef(mod_cause_glm_nb[[i]])$cond+1.96*se)
    axis(1,at=1:ncoef,labels=FALSE)
    text(x=1:ncoef,y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),labels=names(fixef(mod_cause_glm_nb[[i]])$cond),
         srt=45, adj=1, xpd=NA)
}
dev.off()

pdf("graphs/cause_specific_nb_glm_residual_plots_mx.pdf",width=16,height=9)
for (i in 1:length(mod_cause_glm_nb)) {
    X <- model.matrix(mod_cause_glm_nb[[i]]$modelInfo$allForm$formula,data=chn_all)
    preds <- exp(X %*% fixef(mod_cause_glm_nb[[i]])$cond)
    resids <- preds - chn_all$mx
    plot(resids,ylab="Residuals (mx scale)",xlab="Observation",
         main=paste("NegBin GLM:",paste(as.character(mod_cause_glm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
    for (j in dims2) {
        plot(resids~chn_all[,factor(get(j))],xlab=j,ylab="Residuals",
             main=paste("NegBin GLM:",paste(as.character(mod_cause_glm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
    }
}
dev.off()

pdf("graphs/cause_specific_nb_glm_residual_plots_logmx_noInf.pdf",width=16,height=9)
for (i in 1:length(mod_cause_glm_nb)) {
    X <- model.matrix(mod_cause_glm_nb[[i]]$modelInfo$allForm$formula,data=chn_all)
    preds <- X %*% fixef(mod_cause_glm_nb[[i]])$cond
    resids <- preds[!is.infinite(chn_all$logmx)] - chn_all$logmx[!is.infinite(chn_all$logmx)]
    plot(resids,ylab="Residuals (log scale, no -Inf)",xlab="Observation",
         main=paste("NegBin GLM:",paste(as.character(mod_cause_glm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
    for (j in dims2) {
        plot(resids~chn_all[!is.infinite(chn_all$logmx),factor(get(j))],xlab=j,ylab="Residuals",
             main=paste("NegBin GLM:",paste(as.character(mod_cause_glm_nb[[i]]$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
    }
}
dev.off()

## cause correlations in residuals
i <- length(mod_cause_glm_nb)
X <- model.matrix(mod_cause_glm_nb[[i]]$modelInfo$allForm$formula,data=chn_all)
preds <- X %*% fixef(mod_cause_glm_nb[[i]])$cond
chn_all$preds <- preds
chn_all$exp_resids <- exp(chn_all$preds) - chn_all$mx
chn_all$log_resids <- chn_all$preds - chn_all$logmx

cause_resids <- as.data.frame(reshape(chn_all[chn_all$logmx!=-Inf,c("log_resids","exp_resids","cause","agegp_name","year","reg","res")],
                                      direction="wide",v.names=c("log_resids","exp_resids"),timevar=c("cause"),
                                      idvar=c("agegp_name","year","reg","res")))
cor(cause_resids[,grep("exp_resids",names(cause_resids),value=T)],use="pairwise.complete.obs")
cor(cause_resids[,grep("log_resids",names(cause_resids),value=T)],use="pairwise.complete.obs")

pdf("graphs/cause_resid_pairs_plot.pdf",width=16,height = 9)
for (i in (1:length(deathvars))) {
    for (j in 1:i) {
        if (i == j) next
        c1 <- deathvars[i]
        c2 <- deathvars[j]
        plot(cause_resids[,gsub("d_","log_resids.",c1)] ~ cause_resids[,gsub("d_","log_resids.",c2)],
             xlab=c2,ylab=c1)
    }
}
pairs(cause_resids[,grep("log_resids",names(cause_resids),value=T)])
dev.off()

## GLMMs on cause-specific mortality

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
AIC(best.cause.glm.nb,best.cause.glmm.nb)

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

pdf("graphs/cause_resid_pairs_plot_glmm_yearRE_corcauseRE.pdf",width=16,height = 9)
pairs(cause_resids2[,grep("resids",names(cause_resids2),value=T)],
      main=paste("NegBin GLM:",paste(as.character(mod3$modelInfo$allForm$formula)[-2],collapse=" "),sep=" "))
dev.off()

#############
## Testing different dimensions for time RE
#############

## raw data
chn_agg_age_year <- aggregate(cbind(chn_all$deaths,chn_all$exposure),by=list(chn_all$agegp,chn_all$year),sum)
chn_agg_age_year$mx <- chn_agg_age_year[,3]/chn_agg_age_year[,4]
chn_agg_age_year$col <- as.numeric(as.factor(chn_agg_age_year$Group.1))
chn_agg_cause_year <- aggregate(cbind(chn_all$deaths,chn_all$exposure),by=list(chn_all$cause,chn_all$year),sum)
chn_agg_cause_year$mx <- chn_agg_cause_year[,3]/chn_agg_cause_year[,4]
chn_agg_cause_year$col <- as.numeric(as.factor(chn_agg_cause_year$Group.1))
chn_agg_strata_year <- aggregate(cbind(chn_all$deaths,chn_all$exposure),by=list(chn_all$strata,chn_all$year),sum)
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

#############
## Plot model predictions and compare to raw data
#############
causes <- unique(chn_all$cause)
chn_all$preds <- predict(mod7,type="link") - log(chn_all$exposure)
setwd(wd)
pdf("graphs/glmm_model_preds_vs_data_over_time_noagecauseinteraction.pdf",width=16,height=9)
for (i in 1:length(unique(chn_all$strata))) {
    for (j in 1:length(unique(chn_all$cause))) {
        par(mfrow=c(2,3))
        for (k in 1:length(unique(chn_all$agegp))) {
            tmp <- chn_all[chn_all$strata==stratas[i] & chn_all$cause == causes[j] & chn_all$agegp == ages[k], ]
            plot(logmx~year,data=tmp,col=alpha("red",0.5),pch=19,ylim=range(c(tmp$preds,tmp$logmx[!is.infinite(tmp$logmx)])),
                 main=paste0("cause: ",causes[j],"; strata: ",stratas[i],"; age group: ",ages[k]))
            lines(preds~year,data=tmp,col=alpha("blue",0.5),pch=19)
        }
    }
}
dev.off() 
