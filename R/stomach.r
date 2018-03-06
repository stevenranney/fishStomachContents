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
library(scales)
library(ggpubr)

set.seed(256)

source("R/helper_functions.R")

#From Gabelhouse
wae_length_cat <- c(250, 380, 510, 630, 760)
smb_length_cat <- c(180, 280, 350, 430, 510)

#Ws parameters
wae_ws <- c(-5.453, 3.180)
smb_ws <- c(-5.329, 3.2)

#Read in and set up dataframe, keeping only rows related to SMB and WAE
stomach <- 
  read.csv("data/stomach_contents.csv", header=T) %>%
  filter(species %in% c("SMB", "WAE")) %>%
  mutate(weight_empty = weight - st_weight, 
         rel_weight = ifelse(species == "SMB", 
                             calc_rel_weight(weight, length, smb_ws), 
                             calc_rel_weight(weight, length, wae_ws)), 
         rel_weight_empty = ifelse(species == "SMB", 
                                   weight_empty %>% calc_smb_wr(length),
                                   weight_empty %>% calc_wae_wr(length)), 
         psd = ifelse(species == "SMB", 
                      assign_length_cat(length, smb_length_cat), 
                      assign_length_cat(length, wae_length_cat)), 
         psd = psd %>% factor(levels = c("SS", 
                                         "S-Q", 
                                         "Q-P", 
                                         "P-M", 
                                         "M-T", 
                                         ">T")), 
         length_class = length %>% round_down() + 5, 
         lake = lake %>% as.factor())

# Create maximum stomach volume regressions
max_st_contents_models <- 
  stomach %>%
  filter(st_weight > 0) %>%
  group_by(species, lake, psd) %>%
  filter(st_weight == max(st_weight)) %>%
  group_by(species) %>% 
  do(lm = lm(st_weight ~ weight_empty, data = .), 
     rq = rq(st_weight ~ weight_empty, data = ., tau = 0.99), 
     rq_null = rq(st_weight ~ 1, data = ., tau = 0.99)) 

#Write model parameters and goodness of fit values for both models for both species
data.frame(species = c(rep("SMB", 2), rep("WAE", 2)), 
           model = c(rep(c("lm", "rq"), 2)), 
           intercept = c((max_st_contents_models %>% filter(species == "SMB") %>% 
                            pull(lm) %>% pluck(1) %>% coefficients())[[1]], 
                         (max_st_contents_models %>% filter(species == "SMB") %>% 
                              pull(rq) %>% pluck(1) %>% coefficients())[[1]], 
                         (max_st_contents_models %>% filter(species == "WAE") %>% 
                            pull(lm) %>% pluck(1) %>% coefficients())[[1]], 
                         (max_st_contents_models %>% filter(species == "WAE") %>% 
                            pull(rq) %>% pluck(1) %>% coefficients())[[1]]), 
           slope = c((max_st_contents_models %>% filter(species == "SMB") %>% 
                        pull(lm) %>% pluck(1) %>% coefficients())[[2]], 
                     (max_st_contents_models %>% filter(species == "SMB") %>% 
                        pull(rq) %>% pluck(1) %>% coefficients())[[2]], 
                     (max_st_contents_models %>% filter(species == "WAE") %>% 
                        pull(lm) %>% pluck(1) %>% coefficients())[[2]], 
                     (max_st_contents_models %>% filter(species == "WAE") %>% 
                        pull(rq) %>% pluck(1) %>% coefficients())[[2]]),
           gof = c(rep(c("R2", "R1"), 2)), 
           value = c((max_st_contents_models %>% filter(species == "SMB") %>% 
                       pull(lm) %>% pluck(1) %>% summary())$r.squared, 
                      R1(max_st_contents_models %>% filter(species == "SMB") %>% 
                           pull(rq) %>% pluck(1), 
                         max_st_contents_models %>% filter(species == "SMB") %>% 
                           pull(rq_null) %>% pluck(1)), 
                     (max_st_contents_models %>% filter(species == "WAE") %>% 
                        pull(lm) %>% pluck(1) %>% summary())$r.squared, 
                     R1(max_st_contents_models %>% filter(species == "WAE") %>% 
                          pull(rq) %>% pluck(1), 
                        max_st_contents_models %>% filter(species == "WAE") %>% 
                          pull(rq_null) %>% pluck(1)))) %>%
  write.csv("output/model_params_and_gof_values.csv", row.names = FALSE)

# Add estimated max weight stomach contents weight to all individuals
stomach <- 
  stomach %>%
  mutate(weight_max_rq = weight_empty + ifelse(species == "SMB", 
                                               weight_empty %>%
                                                 apply_lm(max_st_contents_models %>% filter(species == "SMB") %>% pull(rq) %>% pluck(1)),
                                               weight_empty %>%
                                                 apply_lm(max_st_contents_models %>% filter(species == "WAE") %>% pull(rq) %>% pluck(1))), 
         rel_weight_max_rq = ifelse(species == "SMB", 
                                    weight_max_rq %>% calc_smb_wr(length),
                                    weight_max_rq %>% calc_wae_wr(length)))


labels <- c("SMB" = "Smallmouth bass", 
            "WAE" = "Walleye", 
            "Clear" = "Clear", 
            "Enemy.Swim" = "Enemy Swim", 
            "MO.river" = "MO River", 
            "Pickerel" = "Pickerel", 
            "Roy" = "Roy", 
            "Bitter" = "Bitter", 
            "Harlan.res" = "Harlan Reservoir", 
            "Kansas" = "Kansas", 
            "Pelican" = "Pelican", 
            "Twin" = "Twin")

stomach %>%
  filter(st_weight > 0) %>%
  group_by(species, lake, psd) %>%
  filter(st_weight == max(st_weight)) %>%
  ggplot(aes(x = weight_empty, y = st_weight)) +
  geom_point() +
  labs(x = expression("Total weight minus stomach contents (" ~ italic(W)[E] ~ "; g)"), 
       y = expression("Stomach contents weight (" ~ italic(W)[St] ~ "; g)")) +
  geom_smooth(aes(linetype = "1"), 
              method = "lm",
              se = FALSE,
              formula = y ~ x,
              colour = "black") +
  geom_smooth(aes(linetype = "2"), 
              method = "rq",
              se = FALSE,
              formula = y ~ x,
              method.args = list(tau = 0.99), 
              colour = "black") +
  scale_x_continuous(labels = comma) +
  facet_wrap(~species, scales = "free", labeller = as_labeller(labels)) +
  scale_linetype_discrete(name = "Model", labels = c("Linear", expression(99^th ~ "Quantile"))) +
  theme_bw() +
  theme(legend.position = "bottom", 
        strip.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("output/model_figure.png")
ggsave("output/model_figure.tiff")

stomach <- 
  stomach %>%
  group_by(species, lake, psd) %>%
  mutate(n = n()) %>%
  filter(n >= 3) # Remove combinations of species x lake x psd where sample size < 3

# Determine if within species within population within psd rel_weight values are normally distributed
stomach %>% 
  group_by(species, lake, psd) %>% 
  summarize(norm_test_pvalue = ifelse(length(rel_weight) > 3, shapiro.test(rel_weight)$p.value, NA)) %>%
  write.csv(paste0("output/", "shaprio_test_rel_weight_by_sp_lake_psd.csv"), row.names = FALSE)

# Determine if within species within population within psd rel_weight values 
# are significantly different from rel_weight_empty or rel_weight_max values
stomach %>% 
  group_by(species, lake, psd) %>% 
  summarize(mw_wr_wre = wilcox.test(rel_weight, y = rel_weight_empty, data = .)$p.value, 
            #mw_wr_wrm_lm = wilcox.test(rel_weight, y = rel_weight_max_lm, data = .)$p.value, 
            mw_wr_wrm_rq = wilcox.test(rel_weight, y = rel_weight_max_rq, data = .)$p.value, 
            #mw_wre_wrm_lm = wilcox.test(rel_weight_empty, y = rel_weight_max_lm, data = .)$p.value, 
            mw_wre_wrm_rq = wilcox.test(rel_weight_empty, y = rel_weight_max_rq, data = .)$p.value) %>%
  write.csv("output/mann_whitney_wr_diffs.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# summary table with Sample size, wilcox.test, and %diff for all species, lake, psd, Wr comb
summary_table <- 
  stomach %>%
  group_by(species, psd, lake) %>%
  summarize(n = n(), 
            wre_wr = wilcox.test(rel_weight_empty, rel_weight)$p.value %>% signif(3), 
            wre_wr_diff = calc_perc_diff(median(rel_weight_empty), median(rel_weight)) %>% signif(3), 
            wre_wrmax = wilcox.test(rel_weight_empty, rel_weight_max_rq)$p.value %>% signif(3), 
            wre_wrmax_diff = calc_perc_diff(median(rel_weight_empty), median(rel_weight_max_rq))%>% signif(3), 
            wr_wrmax = wilcox.test(rel_weight, rel_weight_max_rq)$p.value %>% signif(3), 
            wr_wrmax_diff = calc_perc_diff(median(rel_weight), median(rel_weight_max_rq)) %>% signif(3)) %>%
  mutate(wre_wr = ifelse(wre_wr < 0.05, paste0(wre_wr, "*"), wre_wr %>% as.character()), 
         wre_wrmax = ifelse(wre_wrmax < 0.05, paste0(wre_wrmax, "*"), wre_wrmax %>% as.character()), 
         wr_wrmax = ifelse(wr_wrmax < 0.05, paste0(wr_wrmax, "*"), wr_wrmax %>% as.character()))

summary_table %>%
  write.csv("output/summary_table.csv", row.names = F)

#-------------------------------------------------------------------------------
# Summary table of 25th, 50th, and 75th quantile (first, second, and third quartiles) 
# to go with the plots

wr_values <- 
  stomach %>%
  group_by(species, psd, lake) %>%
  summarize(wr_first = quantile(rel_weight, probs = 0.25), 
            wre_first = quantile(rel_weight_empty, probs = 0.25), 
            wrmax_first = quantile(rel_weight_max_rq, probs = 0.25), 
            wr_second = median(rel_weight), 
            wre_second = median(rel_weight_empty), 
            wrmax_second = median(rel_weight_max_rq),
            wr_third = quantile(rel_weight, probs = 0.75), 
            wre_third = quantile(rel_weight_empty, probs = 0.75), 
            wrmax_third = quantile(rel_weight_max_rq, probs = 0.75))

wr_values[, 1:3] %>%
  bind_cols(round(wr_values[,4:length(wr_values)], 0)) %>%
  write.csv("output/summary_iqrange.csv", row.names = F)


#-------------------------------------------------------------------------------
# Boxplot 

measure_vars <- c("rel_weight", "rel_weight_empty", 
                  "rel_weight_max_rq")#, "rel_weight_max_lm")

stomach %>% 
  filter(psd != ">T") %>%
  melt(id.vars = c("species", "lake", "psd"), 
       measure.vars = measure_vars) %>%
  ggplot(aes(x = psd, y = value, fill = variable)) +
  geom_boxplot(outlier.colour = NA) + #outliers not displayed
  facet_wrap(~species, scales = "free_y", labeller = as_labeller(labels)) +
  coord_cartesian(ylim = c(50, 180)) +
  scale_y_continuous(breaks = seq(60, 180, by = 20)) +
  scale_fill_grey(name = "Relative weight calculation", 
                  labels = c(expression(W[r]), expression(W[rE]), 
                             expression(W[rMax]))) + #, expression(W[rMaxQ]))) +
  theme_bw() +
  theme(legend.position = "bottom", 
        strip.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Length category", y = "Relative weight value")

ggsave("output/boxplot_fig.png")
ggsave("output/boxplot_fig.tiff")

###############################################################################################
# Boxplots by species and population

#SMB
smb_plot <- 
  stomach %>% 
  filter(psd != ">T") %>%
  melt(id.vars = c("species", "lake", "psd"), 
       measure.vars = measure_vars) %>%
  filter(species == "SMB") %>%
  ggplot(aes(x = psd, y = value, fill = variable)) +
  geom_boxplot(outlier.colour = NA) + #outliers not displayed
  facet_wrap(~lake, scales = "free", labeller = as_labeller(labels), 
             drop = TRUE) +
  coord_cartesian(ylim = c(50, 180)) +
  scale_y_continuous(breaks = seq(60, 180, by = 20)) +
  scale_fill_manual(name = "Relative weight calculation", 
                    labels = c(expression(W[r]), expression(W[rE]), 
                               expression(W[rMax])), 
                    values = gray.colors(3, start = 0.4)) + #, expression(W[rMaxQ]))) +
  theme_bw() +
  theme(legend.position = "bottom", 
        plot.margin = margin(0,0.1,0,0.1, "cm"),
        strip.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Length category", y = "Relative weight value")

ggsave("output/smallmouth_wr_plot.png", plot = smb_plot)
ggsave("output/smallmouth_wr_plot.tiff", plot = smb_plot)

#WAE
wae_plot <- 
  stomach %>% 
  filter(psd != ">T") %>%
  melt(id.vars = c("species", "lake", "psd"), 
       measure.vars = measure_vars) %>%
  filter(species == "WAE") %>%
  ggplot(aes(x = psd, y = value, fill = variable)) +
  geom_boxplot(outlier.colour = NA) + #outliers not displayed
  facet_wrap(~lake, scales = "free", labeller = as_labeller(labels), 
             drop = TRUE) +
  # coord_cartesian(ylim = c(50, 180)) +
  scale_fill_manual(name = "Relative weight calculation", 
                  labels = c(expression(W[r]), expression(W[rE]), 
                             expression(W[rMax])), 
                  values = gray.colors(3, start = 0.4)) +# , expression(W[rMaxQ]))) +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.margin = margin(0,0.1,0,0.1, "cm"),
        strip.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Length category", y = "Relative weight value")

ggsave("output/walleye_wr_plot.png", plot = wae_plot)
ggsave("output/walleye_wr_plot.tiff", plot = wae_plot)

# Combine both plots into one
ggarrange(smb_plot + theme(legend.position = "none"), wae_plot,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

ggsave(paste0("output/combine_bw.png"), width = 8.5, height = 11, units = "in")
ggsave(paste0("output/combine_bw.tiff"), width = 8.5, height = 11, units = "in")




