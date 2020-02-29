# Austin Schumacher
# 2/11/2020
# Plot CSMFs in each age-region over time

rm(list=ls())

library(INLA); library(tidyverse); library(gridExtra); library(ggpubr);

# directory where you save your results from '../model_fitting/model_fit_inla.R' (set these yourself)
datadir <- "~/Dropbox/dissertation_2/cause_specific_child_mort/estimation_china/results"
setwd(datadir)

# load data
mod_inla <- readRDS("china_results_pcpriors.RDS")
res <- readRDS("china_results_long.RDS")

# plot cause fractions
ggplot(data = res[res$measure == "pred_csmf_m",], 
       aes(x = year_name, y = value, fill = cause_name)) +
  geom_area(position = 'stack') +
  facet_grid(reg2 ~ agegp_factor) +
  ggtitle(NULL) +
  xlab("year") +
  ylab("cause fraction") +
  labs(fill = "cause") +
  # theme(strip.placement = "outside") +
  theme_light()
ggsave(filename = paste0("../graphs/cause_fractions_pcpriors_all.pdf"), 
       device = "pdf", width = 16, height = 9)

## for publication
ggplot(data = res[res$measure == "pred_csmf_m",], 
       aes(x = year_name, y = value, fill = cause_name)) +
  geom_area(position = 'stack') +
  scale_fill_viridis_d(option = "B") +
  facet_grid(reg2 ~ agegp_factor) +
  ggtitle(NULL) +
  xlab("year") +
  ylab("cause fraction") +
  labs(fill = "") +
  theme_light() + 
  theme(legend.position = "top", 
        legend.text = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 30, vjust = 0.9, hjust = 0.9),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 14, face = "bold", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  guides(fill = guide_legend(label.position = "bottom"))
ggsave(filename = paste0("../graphs/cause_fractions_pcpriors_all_pub.pdf"), 
       device = "pdf", width = 9, height = 9)
