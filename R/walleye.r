
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
  summarise(shapiro_test_pvalue = ifelse(length(rel_weight) < 3, NA, (rel_weight %>% shapiro.test())$p.value)) %>%
  arrange(lake, psd)


#Another method to calculate Coeffecient of determination
#SST <- sum((wae.vol.length$vol-mean(predict(max.vol)))^2)
#SSE <- sum((wae.vol.length$vol-(predict(max.vol)))^2)
#1-(SSE/SST)



#Percent difference between Wr and WrS
((tapply(WrE, psd, median)-tapply(Wr, psd, median))/tapply(Wr, psd, median))*100

#Sample sizes
tapply(wae$rel_weight, wae$psd, FUN = function(x) sum(!is.na(x)))

shapiro.test(wae$rel_weight)

wae %>%
  group_by(psd) %>%
  wilcox.test(rel_weight~rel_weight_empty, data = .)



74with(wae[psd=="1",],
  wilcox.test(Wr, y=WrE))
with(wae[psd=="2",],
  wilcox.test(Wr, y=WrE))
with(wae[psd=="3",],
  wilcox.test(Wr, y=WrE))
with(wae[psd=="4",],
  wilcox.test(Wr, y=WrE))

#Comparing WrE to WrMAX
with(wae[psd=="0",],
    wilcox.test(WrE, y=Wr.max))
with(wae[psd=="1",],
    wilcox.test(WrE, y=Wr.max))
with(wae[psd=="2",],
    wilcox.test(WrE, y=Wr.max))
with(wae[psd=="3",],
  wilcox.test(WrE, y=Wr.max))
with(wae[psd=="4",],
  wilcox.test(WrE, y=Wr.max))

#Percent difference between WrE and Wr.max
((tapply(WrE, psd, median)-tapply(Wr.max, psd, median))/tapply(Wr.max, psd, median))*100

#Comparing Wr to WrMAX
with(wae[psd=="0",],
    wilcox.test(Wr, y=Wr.max))
with(wae[psd=="1",],
    wilcox.test(Wr, y=Wr.max))
with(wae[psd=="2",],
    wilcox.test(Wr, y=Wr.max))
with(wae[psd=="3",],
    wilcox.test(Wr, y=Wr.max))
with(wae[psd=="4",],
    wilcox.test(Wr, y=Wr.max))

#Percent difference between Wr and Wr.max
((tapply(Wr, psd, median)-tapply(Wr.max, psd, median))/tapply(Wr.max, psd, median))*100


###############################################################################################
#Barplot of the whole shebang!
#Adjusted to reflect MEDIANS and quartiles instead of MEAN and 95% CI

z.wae <- matrix(c(tapply(Wr, psd, median),
                tapply(WrE, psd, median),
                tapply(Wr.max, psd, median)),
                nrow=3, ncol=5, byrow=T, dimnames=list(c("Wr", "WrE", "WrMAX"), c("0","1","2","3","4")))

barplot(z.wae, beside=T, space=c(0,1), ylim=c(85,110), xlab="Length Category", axis.lty=1, ylab=expression(paste("Relative weight   ", (italic(W[r])))),
        legend=F, xpd=F, names.arg=c("Substock","S-Q","Q-P","P-M","M-T"))#, density=c(0,0,6), angle=c(0,0,45), col="black")
abline(h=85)
legend(0.5,111, legend=c(expression(italic(W[r]), italic(W[r])[E], italic(W[r])[MAX])), fill=gray.colors(3), bty="n", cex=2)

#Generate matrices of upper and lower CI bars
z.ci.u <- matrix(c(tapply(Wr, psd, median)+((tapply(Wr, psd, sd)/sqrt(tapply(Wr, psd, values)))),
                   tapply(WrE, psd, median)+((tapply(WrE, psd, sd)/sqrt(tapply(WrE, psd, values)))),
                   tapply(Wr.max, psd, median)+((tapply(Wr.max, psd, sd)/sqrt(tapply(Wr.max, psd, values))))),
                   nrow=3, ncol=5, byrow=T, dimnames=list(c("Wr", "WrE", "WrMAX"), c("0","1","2","3","4"))
                 )
z.ci.l <- matrix(c(tapply(Wr, psd, median)-((tapply(Wr, psd, sd)/sqrt(tapply(Wr, psd, values)))),
                   tapply(WrE, psd, median)-((tapply(WrE, psd, sd)/sqrt(tapply(WrE, psd, values)))),
                   tapply(Wr.max, psd, median)-((tapply(Wr.max, psd, sd)/sqrt(tapply(Wr.max, psd, values))))),
                   nrow=3, ncol=5, byrow=T, dimnames=list(c("Wr", "WrE", "WrMAX"), c("0","1","2","3","4"))
                )

#matrices of 25th and 75th percentile
z.q.l.wae <- matrix(c(quantile(Wr[psd==0], 0.25, na.rm=T), quantile(Wr[psd==1], 0.25, na.rm=T), quantile(Wr[psd==2], 0.25, na.rm=T), quantile(Wr[psd==3], 0.25, na.rm=T), quantile(Wr[psd==4], 0.25, na.rm=T),
                    quantile(WrE[psd==0], 0.25, na.rm=T), quantile(WrE[psd==1], 0.25, na.rm=T), quantile(WrE[psd==2], 0.25, na.rm=T), quantile(WrE[psd==3], 0.25, na.rm=T), quantile(WrE[psd==4], 0.25, na.rm=T),
                    quantile(Wr.max[psd==0], 0.25, na.rm=T), quantile(Wr.max[psd==1], 0.25, na.rm=T), quantile(Wr.max[psd==2], 0.25, na.rm=T), quantile(Wr.max[psd==3], 0.25, na.rm=T), quantile(Wr.max[psd==4], 0.25, na.rm=T)),
                    nrow=3, ncol=5, byrow=T, dimnames=list(c("Wr", "WrE", "WrMAX"), c("0","1","2","3","4"))
                   )

z.q.u.wae <- matrix(c(quantile(Wr[psd==0], 0.75, na.rm=T), quantile(Wr[psd==1], 0.75, na.rm=T), quantile(Wr[psd==2], 0.75, na.rm=T), quantile(Wr[psd==3], 0.75, na.rm=T), quantile(Wr[psd==4], 0.75, na.rm=T),
                    quantile(WrE[psd==0], 0.75, na.rm=T), quantile(WrE[psd==1], 0.75, na.rm=T), quantile(WrE[psd==2], 0.75, na.rm=T), quantile(WrE[psd==3], 0.75, na.rm=T), quantile(WrE[psd==4], 0.75, na.rm=T),
                    quantile(Wr.max[psd==0], 0.75, na.rm=T), quantile(Wr.max[psd==1], 0.75, na.rm=T), quantile(Wr.max[psd==2], 0.75, na.rm=T), quantile(Wr.max[psd==3], 0.75, na.rm=T), quantile(Wr.max[psd==4], 0.75, na.rm=T)),
                    nrow=3, ncol=5, byrow=T, dimnames=list(c("Wr", "WrE", "WrMAX"), c("0","1","2","3","4"))
                    )

errbar((c(1.5:3.5, 5.5:7.5, 9.5:11.5, 13.5:15.5, 17.5:19.5)),
        z.wae,
        z.q.u.wae,
        z.q.l.wae,
        add=T, pch=" ")


#detach(wae)
#rm(list=ls(all=TRUE))

#Least-squares regression on log10-transformed data
#wae.vol.length <- read.table("wae.vol.length.txt")
#mod.1 <- lm(log10(wae.vol.length$vol)~log10(wae.vol.length$Length))
#plot(log10(wae.vol.length$vol)~log10(wae.vol.length$Length))
#  abline(mod.1)
#summary(mod.1)

#dffits(mod.1)
#length(wae.vol.length$Length)
#2*(sqrt(3/29))
#wae.vol.length <- cbind(wae.vol.length, abs(dffits(mod.1)))
#length.1 <- ifelse(abs(dffits(mod.1))>0.6432675, NA, wae.vol.length$Length)
#wae.vol.length <- cbind(wae.vol.length, length.1)
#vol.1 <- ifelse(length.1>0, wae.vol.length$vol, NA)
#wae.vol.length <- cbind(wae.vol.length, vol.1)

#mod.2 <- lm(log10(wae.vol.length$vol.1)~log10(wae.vol.length$length.1))
#plot(log10(wae.vol.length$vol.1)~log10(wae.vol.length$length.1))
#  abline(mod.2)
#  abline(mod.1, col="red")
#summary(mod.2)

