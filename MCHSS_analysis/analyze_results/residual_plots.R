# Austin Schumacher
# 2/11/2020
# Plot distributions of residuals from our final model across various combinations of region, age, time, and cause

rm(list=ls())

library(INLA); library(tidyverse); library(gridExtra); library(ggpubr); library(data.table);

# directory for results
savedir <- "../../../../Dropbox/SRS-child-mortality-output/"

# create folders to store results if necessary
if (!file.exists(paste0(savedir, "graphs"))) {
    dir.create(paste0(savedir, "graphs"))
}
if (!file.exists(paste0(savedir, "graphs/residual_plots"))) {
    dir.create(paste0(savedir, "graphs/residual_plots"))
}

# load and format our results
res <- readRDS(paste0(savedir, "results/china_results_long.RDS"))
res <- res %>% 
    mutate(agegp_factor = factor(agegp_name,
                                 levels = c("0-6d", "7-27d", "1-5m", "6-11m", "12-23m", "24-59m")))

## plot hist and density of both types of residuals
pdf(paste0(savedir, "graphs/residual_plots/resid_std_compare_.pdf"),width = 10, height = 6)
ggplot(data = res[res$measure %in% c("resid_std_mu","resid_std_mu_sigma_e"),], 
       aes(x = value, color = measure, fill = measure)) +
  geom_histogram(aes(y=..density..), position = "identity", alpha = 0.3, linetype = "dashed",binwidth = 0.2) +
  geom_density(alpha = 0.2, fill = NA) +
  xlab("residual") +
  theme_light()

plot(res$value[res$measure=="resid_std_mu"] ~ res$value[res$measure=="resid_std_mu_sigma_e"],
     xlab = "residual standardized by muhat*exp(sigma^2_e)",
     ylab = "residual standardized by muhat")
abline(0,1,col="red")
dev.off()

# plotting function for residuals by "facet_var" colored by "color_var"
plot_resid_violin_3_strata <- function(mydata, myresid, by_var, facet_var, strata_var) {
  if (grepl("age",by_var)) {
    leg_title <- "age"
  } else if (grepl("reg",by_var)) {
    leg_title <- "region"
  } else if (grepl("year",by_var)) {
    leg_title <- "year"
  } else if (grepl("cause",by_var)) {
    leg_title <- "cause"
  } else {
    leg_title <- by_var
  }
  mfrow_spec <- n2mfrow(nrow(unique(mydata[mydata$measure == myresid,facet_var])))
  viridis_pal <- "C"
  ggplot(data = mydata[mydata$measure == myresid,], 
         aes(y = value, x = factor(get(by_var)), 
             fill = factor(get(strata_var)))) +
    geom_violin(alpha = 0.6) +
    scale_color_viridis_d(option = viridis_pal) +
    scale_fill_viridis_d(option = viridis_pal) +
    facet_wrap(~ get(facet_var), nrow = mfrow_spec[2], ncol = mfrow_spec[1], scales = "free") +
    labs(fill = leg_title) +
    xlab("year") +
    ylab("residual") +
    theme_light()
}

# plotting function for residuals by "facet_var" colored by "color_var"
plot_resid_violin_2_strata <- function(mydata, myresid, by_var, facet_var) {
  if (grepl("age",by_var)) {
    leg_title <- "age"
  } else if (grepl("reg",by_var)) {
    leg_title <- "region"
  } else if (grepl("year",by_var)) {
    leg_title <- "year"
  } else if (grepl("cause",by_var)) {
    leg_title <- "cause"
  } else {
    leg_title <- by_var
  }
  mfrow_spec <- n2mfrow(nrow(unique(mydata[mydata$measure == myresid,facet_var])))
  viridis_pal <- "C"
  ggplot(data = mydata[mydata$measure == myresid,], 
         aes(y = value, x = factor(get(by_var)), 
             fill = factor(get(by_var)))) +
    geom_violin(alpha = 0.6) +
    scale_color_viridis_d(option = viridis_pal) +
    scale_fill_viridis_d(option = viridis_pal) +
    facet_wrap(~ get(facet_var), nrow = mfrow_spec[2], ncol = mfrow_spec[1], scales = "free") +
    labs(fill = leg_title) +
    xlab("year") +
    ylab("residual") +
    theme(legend.text = element_text(size = 18, face = "bold"),
          plot.title = element_blank(),
          strip.text = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 18, face = "bold"),
          legend.title = element_text(size = 18, face = "bold")) +
    theme_light()
}

# resid over time colored by age
plot_resid_violin_2_strata(mydata = res, 
                             myresid = "resid_std_mu", 
                             by_var = "agegp_factor", 
                             facet_var = "year_name")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_violin_over_year_by_age.pdf"), 
       device = "pdf", width = 16, height = 9)

# resid over time colored by region
plot_resid_violin_2_strata(mydata = res, 
                             myresid = "resid_std_mu", 
                             by_var = "reg2", 
                             facet_var = "year_name")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_violin_over_year_by_reg.pdf"), 
       device = "pdf", width = 16, height = 9)

# resid over age colored by region
plot_resid_violin_2_strata(mydata = res, 
                             myresid = "resid_std_mu", 
                             by_var = "reg2", 
                             facet_var = "agegp_factor")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_violin_over_age_by_reg.pdf"), 
       device = "pdf", width = 16, height = 9)

# resid over time colored by age
plot_resid_violin_2_strata(mydata = res, 
                           myresid = "resid_std_mu", 
                           by_var = "cause_name", 
                           facet_var = "year_name")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_violin_over_year_by_cause.pdf"), 
       device = "pdf", width = 16, height = 9)

# resid over time colored by region
plot_resid_violin_2_strata(mydata = res, 
                           myresid = "resid_std_mu", 
                           by_var = "cause_name", 
                           facet_var = "reg2")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_violin_over_region_by_cause.pdf"), 
       device = "pdf", width = 16, height = 9)

# resid over age colored by region
plot_resid_violin_2_strata(mydata = res, 
                           myresid = "resid_std_mu", 
                           by_var = "cause_name", 
                           facet_var = "agegp_factor")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_violin_over_age_by_cause.pdf"), 
       device = "pdf", width = 16, height = 9)

# resid over time by age colored by reg
plot_resid_violin_3_strata(mydata = res, 
                             myresid = "resid_std_mu", 
                             by_var = "agegp_factor", 
                             facet_var = "year_name",
                             strata_var = "reg2")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_violin_over_year_by_age_strata_reg.pdf"), 
       device = "pdf", width = 30, height = 10)

# resid over age colored by time
plot_resid_violin_3_strata(mydata = res, 
                             myresid = "resid_std_mu", 
                             by_var = "agegp_factor", 
                             facet_var = "year_name",
                             strata_var = "cause_name")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_violin_over_year_by_age_strata_cause.pdf"), 
       device = "pdf", width = 30, height = 10)

# resid over age colored by time
plot_resid_violin_3_strata(mydata = res, 
                           myresid = "resid_std_mu", 
                           by_var = "reg2", 
                           facet_var = "year_name",
                           strata_var = "cause_name")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_violin_over_year_by_reg_strata_cause.pdf"), 
       device = "pdf", width = 30, height = 10)

# plotting function for residuals by "facet_var" colored by "color_var"
plot_resid_time_2_strata <- function(mydata, myresid, by_var, facet_var) {
  if (grepl("age",by_var)) {
    leg_title <- "age"
  } else if (grepl("reg",by_var)) {
    leg_title <- "region"
  } else if (grepl("year",by_var)) {
    leg_title <- "year"
  } else if (grepl("cause",by_var)) {
    leg_title <- "cause"
  } else {
    leg_title <- by_var
  }
  mfrow_spec <- n2mfrow(nrow(unique(mydata[mydata$measure == myresid,facet_var])))
  viridis_pal <- "C"
  ggplot(data = mydata[mydata$measure == myresid,], 
         aes(y = value, x = year_name, 
             color = factor(get(by_var)))) +
    geom_point(alpha = 0.6) +
    geom_smooth(se = FALSE) +
    scale_color_viridis_d(option = viridis_pal) +
    scale_fill_viridis_d(option = viridis_pal) +
    facet_wrap(~ get(facet_var), nrow = mfrow_spec[2], ncol = mfrow_spec[1], scales = "free") +
    labs(color = leg_title) +
    xlab("year") +
    ylab("residual") +
    theme_light()
}

plot_resid_time_2_strata(mydata = res, 
                           myresid = "resid_std_mu", 
                           by_var = "reg2", 
                           facet_var = "agegp_factor")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_time_over_age_col_reg.pdf"), 
       device = "pdf", width = 16, height = 9)

plot_resid_time_2_strata(mydata = res, 
                         myresid = "resid_std_mu", 
                         by_var = "cause_name", 
                         facet_var = "agegp_factor")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_time_over_age_col_cause.pdf"), 
       device = "pdf", width = 16, height = 9)

plot_resid_time_2_strata(mydata = res, 
                         myresid = "resid_std_mu", 
                         by_var = "cause_name", 
                         facet_var = "reg2")
ggsave(filename = paste0(savedir, "graphs/residual_plots/resid_std_time_over_region_col_cause.pdf"), 
       device = "pdf", width = 16, height = 9)
