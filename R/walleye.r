
library(dplyr)
library(ggplot2)
library(quantreg)

source("R/helper_functions.R")

stomach <- 
  read.csv("data/stomach_contents.csv", header=T) %>%
  filter(species %in% c("WAE", "SMB"))

################################################################################
#For Walleye data
wae <- 
  stomach %>%
  filter(species=="WAE")

wae <- 
  wae %>%
  mutate(rel_weight = calc_wae_wr(weight, length), 
         rel_weight_empty = (weight - st_weight) %>% calc_wae_wr(length), 
         psd = ifelse((length>=250)&(length<380), "S-Q",
                      ifelse((length>=380)&(length<510), "Q-P",
                             ifelse((length>=510)&(length<630), "P-M",
                                    ifelse((length>=630)&(length<760), "M-T",
                                           ifelse(length>=760, ">T", "substock"))))), 
         psd = psd %>% factor(levels = c("substock", 
                                         "S-Q", 
                                         "Q-P", 
                                         "P-M", 
                                         "M-T", 
                                         ">T")), 
         length_class = length %>% round_down() + 5, 
         lake = lake %>% as.factor(), 
         vol = 1.05 * st_weight)

max_vol <- 
  wae %>% 
  group_by(lake, psd) %>%
  filter(vol == max(vol)) %>%
  filter(vol > 0)

max_vol_nls <- 
  max_vol %>% 
  nls(vol~a*length^b, data = ., start=list(a=0.00000005, b=3.0), trace=TRUE,
      control=nls.control(maxiter=10000))

max_vol %>%
  ungroup() %>%
  filter(vol > 0) %>%
  ggplot(aes(x = length, y = vol)) +
  geom_point() + 
  geom_smooth(method = "nls", 
              formula = y ~ a * x^b, 
              se = FALSE, 
              method.args = list(start = list(a = 0.00000005, b = 3.0), 
                                 trace = TRUE, 
                                 control = nls.control(maxiter = 10000)), 
              colour = "black") +
  labs(x = "Total length (mm)", y = "Maximum stomach contents (ml)") +
  ggtitle("nls of maximum stomach volume by length category")

ggsave(paste0("output/", Sys.Date(), "_max_vol_regression.png"))

#Sink normality tests to a log file
sink(paste0("output/", Sys.Date(), "_output.log"))
print(paste0("R^2 for non-linear model of max stomach contents~total length of WAE = ", 
             R2(max_vol$length, max_vol$vol, max_vol_nls)))
sink(type = "message")
sink()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

wae <- 
  wae %>%
  mutate(weight_empty = weight - st_weight)


mod_95 <- 
  wae %>% 
  nlrq(vol~a*length^b, data = ., start = list(a = .00000001, b = 4), tau = .95, 
       trace = TRUE)

wae <- 
  wae %>%
  mutate(mod_95_predict = predict(mod_95))

R2(wae$length, wae$vol, mod_95)

wae %>%
  ungroup() %>%
  filter(vol > 0) %>%
  ggplot(aes(x = length, y = vol)) +
  geom_point() + 
  geom_line(aes(x = length, y = mod_95_predict)) +
  labs(x = "Total length (mm)", y = "Stomach contents (ml)") +
  ggtitle("95th Quantile regresion")

ggsave(paste0("output/", Sys.Date(), "_95th_vol_nlrq_regression.png"))


#Estimate maximum stomach contents weight
wae <- 
  wae %>%
  mutate(max_st_weight_rq = (coef(mod_95)[1]*wae$length^coef(mod_95)[2])/1.05, 
         max_st_weight_nls = (coef(max_vol_nls)[1]*wae$length^coef(max_vol_nls)[2])/1.05,
         rel_weight_max_rq = calc_wae_wr((wae$weight_empty + wae$max_st_weight_rq), wae$length), 
         rel_weight_max_nls = calc_wae_wr((wae$weight_empty + wae$max_st_weight_nls), wae$length))

sink(paste0("output/", Sys.Date(), "_output.log"), append = TRUE)
#t-test for significant differences between Wr, WrE, and WrMax
#NOT FINISHED
t.test(wae$rel_weight, wae$rel_weight_empty)
t.test(wae$rel_weight_empty, wae$WrMax)
t.test(wae$rel_weight_empty, wae$WrMax)
t.test(wae$rel_weight, wae$rel_weight_max_rq)
t.test(wae$rel_weight, wae$rel_weight_max_nls)
sink(type = "message")
sink()



wae %>%
  group_by(psd) %>%
  lm(rel_weight~lake, data = .) %>% summary()

sink(paste0("output/", Sys.Date(), "_output.log"), append = TRUE)
print("substock")
Wr0 <- lm(wae$rel_weight[wae$psd == "substock"]~wae$lake[wae$psd == "substock"])
plot(wae$rel_weight[wae$psd == "substock"]~wae$lake[wae$psd == "substock"])
summary(Wr0)
TukeyHSD(aov(Wr0))

print("S-Q")
Wr1 <- lm(wae$rel_weight[wae$psd == "S-Q"]~wae$lake[wae$psd == "S-Q"])
plot(wae$rel_weight[wae$psd == "S-Q"]~wae$lake[wae$psd == "S-Q"])
summary(Wr1)
TukeyHSD(aov(Wr1))

print("Q-P")
Wr2 <- lm(wae$rel_weight[wae$psd == "Q-P"]~wae$lake[wae$psd == "Q-P"])
plot(wae$rel_weight[wae$psd == "Q-P"]~wae$lake[wae$psd == "Q-P"])
summary(Wr2)
TukeyHSD(aov(Wr2))

print("P-M")
Wr3 <- lm(wae$rel_weight[wae$psd == "P-M"]~wae$lake[wae$psd == "P-M"])
plot(wae$rel_weight[wae$psd == "P-M"]~wae$lake[wae$psd == "P-M"])
summary(Wr3)
TukeyHSD(aov(Wr3))

print("M-T")
Wr4 <- lm(wae$rel_weight[wae$psd == "M-T"]~wae$lake[wae$psd == "M-T"])
plot(wae$rel_weight[wae$psd == "M-T"]~wae$lake[wae$psd == "M-T"])
summary(Wr4)
TukeyHSD(aov(Wr4))
sink(type = "message")
sink()


#Test for normality within length category within population for all combinations
normality_tests <- 
  wae %>%
  group_by(psd, lake) %>%
  summarise(shapiro_test_pvalue = ifelse(length(rel_weight) < 3, NA, (rel_weight %>% shapiro.test())$p.value), 
            sample_sizes = sum(!is.na(rel_weight))) %>%
  arrange(lake, psd)

wae %>%
  group_by(psd) %>%
  summarise(mw_wr_wre = (wilcox.test(rel_weight, y = rel_weight_empty)$p.value), 
            mw_wre_wrm = (wilcox.test(rel_weight_empty, y = rel_weight_max_nls)$p.value), 
            mw_wr_wrm = (wilcox.test(rel_weight, y = rel_weight_max_nls)$p.value))

#Percent difference between WrE and Wr.max
((tapply(WrE, psd, median)-tapply(Wr.max, psd, median))/tapply(Wr.max, psd, median))*100

#Percent difference between Wr and Wr.max
((tapply(Wr, psd, median)-tapply(Wr.max, psd, median))/tapply(Wr.max, psd, median))*100


###############################################################################################
#Barplot of the whole shebang!
#Adjusted to reflect MEDIANS and quartiles instead of MEAN and 95% CI

err_bars <- 
  wae %>%
  filter(psd != ">T") %>%
  melt(id.vars = c("species", "lake", "psd"), 
       measure.vars = c("rel_weight", "rel_weight_empty", 
                        "rel_weight_max_rq", "rel_weight_max_nls")) %>%
  group_by(psd, variable) %>%
  mutate(upper_quart = quantile(value, 0.75, type = 6), 
         lower_quart = quantile(value, 0.25, type = 6), 
         upper_conf = mean(value) + (1.96 * (sd(value))), 
         lower_conf = mean(value) - (1.96 * (sd(value))))


wae %>% 
  filter(psd != ">T") %>%
  melt(id.vars = c("species", "lake", "psd"), 
       measure.vars = c("rel_weight", "rel_weight_empty", 
                        "rel_weight_max_rq", "rel_weight_max_nls")) %>%
  ggplot(aes(x = psd, y = value, fill = variable)) +
  stat_summary(fun.y = "median", geom = "bar", position = "dodge") +
  geom_errorbar(data = err_bars,
                aes(x = psd, ymin = lower_quart, ymax = upper_quart),
                position = "dodge") +
  coord_cartesian(ylim = c(80, 110)) +
  scale_fill_grey(name = "Relative weight value", 
                  #values = gray.colors(5) %>% rev, 
                  labels = c(expression(W[r]), expression(W[rE]), 
                             expression(W[rMaxQ]), expression(W[rMaxNLS]))) +
  theme_bw() +
  theme(legend.position = "bottom", 
        panel.border = element_blank(), 
        axis.line = element_line(color = 'black')) +
  labs(x = "Length category", y = "Relative weight value")
