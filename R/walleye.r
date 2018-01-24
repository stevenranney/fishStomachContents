
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
  labs(x = "Total length (mm)", y = "Maximum stomach contents (ml)")

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
  labs(x = "Total length (mm)", y = "Maximum stomach contents (ml; from quantile regression)")


#Estimate maximum stomach contents weight
wae$maxStWeight <- coef(quantMod95)[1]*wae$deltaWeight^coef(quantMod95)[2]

wae$WrMax <- (((wae$Weight-wae$St.weight)+wae$maxStWeight)/(10^(-5.453+3.18*(log10(wae$Length)))))*100

#t-test for significant differences between Wr, WrE, and WrMax
#NOT FINISHED
t.test(wae$Wr, wae$WrE)
t.test(wae$WrE, wae$WrMax)
t.test(wae$Wr, wae$WrMax)

Wr0 <- lm(wae$Wr[wae$psd == "0"]~wae$Lake[wae$psd == "0"])
plot(wae$Wr[wae$psd == "0"]~wae$Lake[wae$psd == "0"])
summary(Wr0)
TukeyHSD(aov(Wr0))

Wr1 <- lm(wae$Wr[wae$psd == "1"]~wae$Lake[wae$psd == "1"])
plot(wae$Wr[wae$psd == "1"]~wae$Lake[wae$psd == "1"])
summary(Wr1)
TukeyHSD(aov(Wr1))

Wr2 <- lm(wae$Wr[wae$psd == "2"]~wae$Lake[wae$psd == "2"])
plot(wae$Wr[wae$psd == "2"]~wae$Lake[wae$psd == "2"])
summary(Wr2)
TukeyHSD(aov(Wr2))

Wr3 <- lm(wae$Wr[wae$psd == "3"]~wae$Lake[wae$psd == "3"])
plot(wae$Wr[wae$psd == "3"]~wae$Lake[wae$psd == "3"])
summary(Wr3)
TukeyHSD(aov(Wr3))

Wr4 <- lm(wae$Wr[wae$psd == "4"]~wae$Lake[wae$psd == "4"])
plot(wae$Wr[wae$psd == "4"]~wae$Lake[wae$psd == "4"])
summary(Wr4)
TukeyHSD(aov(Wr4))


#Test for normality within length category within population for all combinations

normTest <- NULL
for(i in 1:length(unique(wae$psd))){
  for(j in 1:length(unique(wae$Lake))){
      if(values(wae$Wr[wae$psd == unique(wae$psd)[i] & wae$Lake == unique(wae$Lake)[j]]) > 3){
        tmp <- shapiro.test(wae$Wr[wae$psd == unique(wae$psd)[i] & wae$Lake == unique(wae$Lake)[j]])
        normTest <- c(normTest, cbind(as.character(unique(wae$psd))[i], 
                                      as.character(unique(wae$Lake))[j], 
                                      as.numeric(tmp$p.value)))
      } else {
        normTest <- c(normTest, cbind(as.character(unique(wae$psd))[i], 
                                      as.character(unique(wae$Lake))[j], 
                                      "NA"))
      } 
  }
}

normTest <- as.data.frame(matrix(normTest, ncol=3, nrow=30, byrow = TRUE))
#  aggregate(wae$Wr, by=list(wae$psd, wae$Lake), values)
normTest <- normTest[order(normTest[2], normTest[1]), ]

#Another method to calculate Coeffecient of determination
#SST <- sum((wae.vol.length$vol-mean(predict(max.vol)))^2)
#SSE <- sum((wae.vol.length$vol-(predict(max.vol)))^2)
#1-(SSE/SST)


#Calculate maximum stomach weight of individual fish
max.st.weight <- (2.200e-10*Length^4.121)/1.05
wae <- cbind(wae, max.st.weight)

#Calculate Wr.max of individual fish (i.e., total weight minus stomach weight
#plus estimated maximum stomach weight)
Wr.max <- ((Weight-St.weight+max.st.weight)/(10^(-5.453+3.18*(log10(Length))))*100)
wae <- cbind(wae, Wr.max)

#Percent difference between Wr and WrS
((tapply(WrE, psd, median)-tapply(Wr, psd, median))/tapply(Wr, psd, median))*100

#Function to determine sample size
values=function(x)
  {
  return(sum(!is.na(x)))
  }
tapply(Wr, psd, values)

shapiro.test(Wr)

with(wae[psd=="0",],
  wilcox.test(Wr, y=WrE))
with(wae[psd=="1",],
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

