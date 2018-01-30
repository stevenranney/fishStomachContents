################################################################################
# Code to accompany the manuscript 
# The Effect of Stomach Contents on the Relative Weight (Wr) 
# of Smallmouth Bass and Walleye
# S. H. Ranney, A. Zale, and J. Syslo, 2018
################################################################################ 

library(dplyr)
library(ggplot2)
library(quantreg)
library(reshape2)
library(broom)
library(tidyr)
library(purrr)

set.seed(256)

source("R/helper_functions.R")

stomach <- 
  read.csv("data/stomach_contents.csv", header=T) %>%
  filter(species %in% c("SMB", "WAE")) %>%
  mutate(weight_empty = weight - st_weight, 
         rel_weight = ifelse(species == "SMB", 
                             calc_smb_wr(weight, length), 
                             calc_wae_wr(weight, length)), 
         rel_weight_empty = ifelse(species == "SMB", 
                                   weight_empty %>% calc_smb_wr(length),
                                   weight_empty %>% calc_wae_wr(length)), 
         psd = ifelse(species == "SMB", 
                      assign_smb_psd(length), 
                      assign_wae_psd(length)), 
         psd = psd %>% factor(levels = c("Substock", 
                                         "S-Q", 
                                         "Q-P", 
                                         "P-M", 
                                         "M-T", 
                                         ">T")), 
         length_class = length %>% round_down() + 5, 
         lake = lake %>% as.factor())

max_st_contents_models <- 
  stomach %>%
  filter(st_weight > 0) %>%
  group_by(species, lake, psd) %>%
  filter(st_weight == max(st_weight)) %>%
  group_by(species) %>% 
  do(lm = lm(st_weight ~ weight_empty, data = .), 
     rq = rq(st_weight ~ weight_empty, data = ., tau = 0.95), 
     rq_null = rq(st_weight ~ 1, data = ., tau = 0.95)) 


# Add estimated max weight stomach contents weight to all individuals
stomach <- 
  stomach %>%
  mutate(weight_max_lm = weight_empty + ifelse(species == "SMB", 
                                               weight_empty %>%
                                                 apply_lm(max_st_contents_models %>% filter(species == "SMB") %>% pull(lm) %>% pluck(1)),
                                               weight_empty %>%
                                                 apply_lm(max_st_contents_models %>% filter(species == "WAE") %>% pull(lm) %>% pluck(1))), 
         weight_max_rq = weight_empty + ifelse(species == "SMB", 
                                               weight_empty %>%
                                                 apply_lm(max_st_contents_models %>% filter(species == "SMB") %>% pull(rq) %>% pluck(1)),
                                               weight_empty %>%
                                                 apply_lm(max_st_contents_models %>% filter(species == "WAE") %>% pull(lm) %>% pluck(1))), 
         rel_weight_max_lm = ifelse(species == "SMB", 
                                    weight_max_lm %>% calc_smb_wr(length),
                                    weight_max_lm %>% calc_wae_wr(length)), 
         rel_weight_max_rq = ifelse(species == "SMB", 
                                    weight_max_rq %>% calc_smb_wr(length),
                                    weight_max_rq %>% calc_wae_wr(length)))


labels <- c("SMB" = "Smallmouth bass", 
            "WAE" = "Walleye")

stomach %>%
  filter(st_weight > 0) %>%
  group_by(species, lake, psd) %>%
  filter(st_weight == max(st_weight)) %>%
  ggplot(aes(x = weight_empty, y = st_weight)) +
  geom_point() +
  labs(x = "Weight - stomach contents (g)", y = "Stomach contents weight (g)") +
  geom_smooth(aes(linetype = "1"), 
              method = "lm",
              se = FALSE,
              formula = y ~ x,
              colour = "black") +
  geom_smooth(aes(linetype = "2"), 
              method = "rq",
              se = FALSE,
              formula = y ~ x,
              method.args = list(tau = 0.95), 
              colour = "black") +
  facet_wrap(~species, scales = "free", labeller = as_labeller(labels)) +
  scale_linetype_discrete(name = "Model", labels = c("Linear", expression(95^th ~ "Quantile"))) +
  theme_bw() +
  theme(legend.position = "bottom", 
        strip.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste0("output/", "model_figure.png"))

#Evaluate if within species within population within psd rel_weight values are normally distributed
stomach %>% 
  group_by(species, lake, psd) %>% 
  summarize(norm_test_pvalue = ifelse(length(rel_weight) > 3, shapiro.test(rel_weight)$p.value, NA)) %>%
  write.csv(paste0("output/", "sp_lake_psd_norm_test.csv"), row.names = FALSE)

#Evaluate if within species within population within psd rel_weight values 
# are significantly different from rel_weight_empty or rel_weight_max values
stomach %>% 
  group_by(species, lake, psd) %>% 
  summarize(mw_wr_wre = wilcox.test(rel_weight, y = rel_weight_empty, data = .)$p.value, 
            mw_wr_wrm_lm = wilcox.test(rel_weight, y = rel_weight_max_lm, data = .)$p.value, 
            mw_wr_wrm_rq = wilcox.test(rel_weight, y = rel_weight_max_rq, data = .)$p.value, 
            mw_wre_wrm_lm = wilcox.test(rel_weight_empty, y = rel_weight_max_lm, data = .)$p.value, 
            mw_wre_wrm_rq = wilcox.test(rel_weight_empty, y = rel_weight_max_rq, data = .)$p.value) %>%
  write.csv(paste0("output/", "mann_whitney_wr_diffs.csv"), row.names = FALSE)



###############################################################################################
#Barplot 
#Adjusted to reflect MEDIANS and quartiles instead of MEAN and 95% CI

measure_vars <- c("rel_weight", "rel_weight_empty", 
                  "rel_weight_max_lm", "rel_weight_max_rq")

err_bars <- 
  stomach %>%
  filter(psd != ">T") %>%
  melt(id.vars = c("species", "lake", "psd"), 
       measure.vars = measure_vars) %>%
  group_by(species, psd, variable) %>%
  mutate(upper_quart = quantile(value, 0.75, type = 6), 
         lower_quart = quantile(value, 0.25, type = 6), 
         upper_conf = mean(value) + (1.96 * (sd(value))), 
         lower_conf = mean(value) - (1.96 * (sd(value))))


stomach %>% 
  filter(psd != ">T") %>%
  melt(id.vars = c("species", "lake", "psd"), 
       measure.vars = measure_vars) %>%
  ggplot(aes(x = psd, y = value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) + 
  # stat_summary(fun.y = "median", geom = "boxplot", position = "dodge") +
  # geom_errorbar(data = err_bars,
  #               aes(x = psd, ymin = lower_quart, ymax = upper_quart),
  #               position = "dodge") +
  facet_wrap(~species, scales = "free_y", labeller = as_labeller(labels)) +
  coord_cartesian(ylim = c(50, 180)) +
  scale_fill_grey(name = "Relative weight calculation", 
                  #values = gray.colors(5) %>% rev, 
                  labels = c(expression(W[r]), expression(W[rE]), 
                             expression(W[rMax]), expression(W[rMaxQ]))) +
  theme_bw() +
  theme(legend.position = "bottom", 
        strip.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Length category", y = "Relative weight value")

ggsave(paste0("output/", Sys.Date(), "_walleye_fig.png"))
