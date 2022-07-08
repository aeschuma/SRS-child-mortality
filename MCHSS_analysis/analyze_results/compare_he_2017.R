# Austin Schumacher
# 2/14/2020
# Compare our CSMF estimates to He et al (2017) and the observed CSMFs

# preamble ####
rm(list=ls())

library(INLA); library(tidyverse); library(gridExtra); library(data.table);
library(ggpubr); library(haven); library(cowplot); library(haven);

# directory for results
savedir <- "../../../../Dropbox/SRS-child-mortality-output/"

# create folders to store results if necessary
if (!file.exists(paste0(savedir, "graphs"))) {
    dir.create(paste0(savedir, "graphs"))
}
if (!file.exists(paste0(savedir, "graphs/he_compare"))) {
    dir.create(paste0(savedir, "graphs/he_compare"))
}

# load our results
mod_inla <- readRDS(paste0(savedir, "results/china_results_pcpriors.RDS"))
res <- readRDS(paste0(savedir, "results/china_results_long.RDS"))

# load He et al results (email Li Liu lliu26@jhu.edu for access)
he_res_long <- readRDS("../../data/he_2017_test_data.RDS")

# causes that are shared between He and us
common_causes <- c("prematurity", 
                   "birth asphyxia/trauma", 
                   "injuries", 
                   "diarrhea", 
                   "congenital anomalies",
                   "acute resp. infections")

# aggregate our results to nn and pnn so we can compare death counts
res_agg <- res %>%
  select(year_name, agegp, reg2, cause_name, measure, value) %>%
  mutate(year = year_name) %>%
  mutate(cause_name = as.character(cause_name)) %>%
  filter(measure %in% c("fe_rw_m","deaths","exposure")) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(death_est = exp(fe_rw_m) * exposure,
         age_agg = ifelse(agegp %in% c(1,2), "nn", "pn")) %>%
  group_by(year, reg2, age_agg, cause_name) %>%
  summarise(death_est_agg = sum(death_est),
            deaths_agg = sum(deaths)) %>%
  ungroup() %>%
  group_by(year, age_agg, reg2) %>%
  mutate(csmf_est = death_est_agg/sum(death_est_agg),
         csmf = deaths_agg/sum(deaths_agg)) %>%
  filter(cause_name %in% common_causes) %>%
  select(year, age_agg, reg2, cause_name, csmf, csmf_est)

# for He 2017 comparison data, select only causes in common 
he_res_long <- he_res_long %>% filter(cause_name %in% common_causes)

allres <- full_join(res_agg, he_res_long, by = c("year","age_agg","reg2","cause_name"))

# plot comparison function
plot_csmfs <- function(mydata, xvar, yvar, myxlab, myylab) {
  ggplot(data = mydata,
         aes(x = get(xvar), y = get(yvar), color = cause_name)) +
    geom_point(alpha = 0.4, cex = 4) +
    geom_abline(intercept = 0, slope = 1, 
                color = "black", linetype = "dashed", size = 1.25) +
    scale_color_viridis_d(option = "D") +
    xlab(myxlab) +
    ylab(myylab) +
    labs(color = "cause") +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme_light() +
    theme(legend.text = element_text(size = 18),
          plot.title = element_blank(),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 16, face = "bold"),
          legend.title = element_text(size = 18, face = "bold"))
}

# plot our results compared  to He et al
p1 <- plot_csmfs(mydata = allres, yvar = "csmf_est", xvar = "csmf_he",
                 myylab = "our estimates", myxlab = "He et al (2017)")

# plot our results compared to empirical
p2 <- plot_csmfs(mydata = allres, yvar = "csmf_est", xvar = "csmf",
                 myylab = "our estimates", myxlab = "observed data")

# plot He results compared to empirical
p3 <- plot_csmfs(mydata = allres, yvar = "csmf_he", xvar = "csmf",
                 myxlab = "observed data", myylab = "He et al (2017)")

# save the  plots
myplots <- plot_grid(p1 + theme(legend.position="none"), 
                     p2 + theme(legend.position="none"), 
                     p3 + theme(legend.position="none"),
                     labels = "AUTO", nrow = 1)
mylegend <- get_legend(p1 + theme(legend.position="bottom",
                                  legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(myplots, mylegend, nrow = 2, rel_heights = c(3,0.5))
ggsave(filename = paste0(savedir, "graphs/he_compare/he_csmf_compare.pdf"), 
       device = "pdf",width = 13, height = 4)

# plot csmfs over time function
plot_csmfs_over_time <- function(mydata, mymeasure) {
  ggplot(data = mydata, 
         aes(x = year, y = get(mymeasure), fill = cause_name)) +
    geom_area(position = 'stack') +
    facet_grid(reg2 ~ age_agg) +
    ggtitle(NULL) +
    xlab("year") +
    ylab("cause fraction") +
    labs(fill = "cause") +
    theme_light() +
    guides(fill = guide_legend(label.position = "bottom", nrow = 1)) +
    theme(legend.text = element_text(size = 16),
          plot.title = element_blank(),
          axis.text.x = element_text(size = 16, angle = 30, hjust = 0.9, vjust = 0.9),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 16, face = "bold"),
          legend.title = element_text(size = 16, face = "bold"))
}

# plot estimates
g_csmf <- plot_csmfs_over_time(mydata = allres, mymeasure = "csmf")
g_csmf_est <- plot_csmfs_over_time(mydata = allres, mymeasure = "csmf_est")
g_csmf_he <- plot_csmfs_over_time(mydata = allres, mymeasure = "csmf_he")

prow <- plot_grid(g_csmf + theme(legend.position="none"), 
                  g_csmf_est + theme(legend.position="none"), 
                  g_csmf_he+ theme(legend.position="none"),
                  labels = "AUTO", nrow = 1)
legend <- get_legend(g_csmf + theme(legend.box.margin = margin(0, 0, 0, 12),
                                    legend.position="bottom"))
plot_grid(prow, legend, nrow = 2, rel_heights = c(3,0.5))
ggsave(filename = paste0(savedir, "graphs/he_compare/cause_fractions_emp_est_he_compare.pdf"), 
       device = "pdf",width = 12, height = 7)
