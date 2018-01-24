

library(dplyr)
library(ggplot2)

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
         length_class = length %>% round_down() + 5)

wae %>%
  ggplot(aes(x = rel_weight)) +
  geom_histogram(bins = 80, fill = "red", alpha = 0.5) +
  geom_histogram(aes(x = rel_weight_empty), 
                 bins = 80, fill = "blue", alpha = 0.5)


#Sink normality tests to a log file
sink("output/normality_tests_wae.log")
shapiro.test(wae$rel_weight)

wae %>% 
  filter(psd == "substock") %>%
  kruskal.test(rel_weight~rel_weight_empty, data = .)
wae %>%
  filter(psd=="S-Q") %>%
  kruskal.test(rel_weight~rel_weight_empty, data = .)
wae %>%
  filter(psd=="Q-P") %>%
  kruskal.test(rel_weight~rel_weight_empty, data = .)
wae %>% 
  filter(psd=="P-M") %>%
  kruskal.test(rel_weight~rel_weight_empty, data = .)
wae %>%
  filter(psd == "M-T") %>%
  kruskal.test(rel_weight~rel_weight_empty, data = .)

sink(type = "message")
sink()


###Max Volume based on 10 mm length class
#cls.1 <- as.numeric(as.character(cut(length,seq(120,710,by=9.99),labels=seq(125,705,by=10),include.lowest=TRUE)))
#wae <- cbind(wae, cls.1)
#tapply(vol, cls.1, max)
#x.1 <- c(125,135,145,155,165,175,185,195,205,215,225,235,245,255,265,275,295,305,
#  315,325,335,345,355,365,375,385,395,405,415,425,435,445,455,465,475,485,495,505,
#  515,525,535,545,555,565,575,585,645,655,665,675,695,705)
#y.1 <- c(0.242865,0.646065,1.087170,1.160145,0.005880,1.373400,1.647765,1.790775,
#  1.184925,0.739620,2.056950,6.394500,0.748650,2.215500,2.877000,0.703500,2.545200,18.893700,
#  4.762065,6.116355,3.780000,4.833150,25.473000,8.313270,18.732000,8.169000,13.690635,
#  14.007000,16.568160,33.064500,21.453705,20.496000,11.847150,27.594000,11.534460,23.953650, 
#  12.474000,22.123500,8.991360,16.264500,17.713500,3.658725,24.209325,0.661500,8.676885,
#  1.605555,84.530250,156.310350,7.549500,8.344140,151.713660,136.920000)
#plot(y.1~x.1, pch=19, ylab="Maximum stomach contents (ml)", xlab="length class (10 mm)", 
#  ylim=c(0,200), xlim=c(100, 800),yaxs="i", xaxs="i", las=1, bty="n")
#mod.1 <- nls(y.1~a*x.1^b, start=list(a=0.000000000003, b=4.7), trace=T, control=nls.control(maxiter=1000))
#  mod <- seq(min(length), 800, 0.01)
#  lines(mod, predict(mod.1, list(x.1 = mod)))   
#  text(350, 150, expression(V==2.69e-21*(length)^7.96))
#  text(350, 140, expression(R^2==0.56))
#summary(mod.1)
#  1-(deviance(mod.1)/((length(x.1)-1)*var(y.1)))

vol <- 1.05*st_weight
wae <- cbind(wae, vol)

tapply(vol, cls.5, max)

x <- c(125,175,225,275,325,375,425,475,525,575,625,675,725)
y <- c(1.087170,1.790775,6.394500,2.877000,18.893700,25.473000,33.064500,27.594000,22.123500,24.209325,84.530250,156.310350,136.920000)
plot(y~x, pch=19, ylab="Maximum stomach contents (ml)", xlab="length class (50 mm)", 
     ylim=c(0,200), xlim=c(100, 800),yaxs="i", xaxs="i", las=1, bty="n")
max.vol <- nls(y~a*x^b, start=list(a=0.00000001, b=3.35), trace=T, )
mod <- seq(min(length), 800, 0.01)
lines(mod, predict(max.vol, list(x = mod)))   
text(350, 150, expression(V==1.464e-10*(length)^4.20))
text(350, 140, expression(R^2==0.87))
summary(max.vol)

#R2 value for the model
#Defines the total sum of squares for both models (totalSS=(n-1)*var(y))
totalSS <- (length(x)-1)*var(y)
#residSS=deviance(nls model)
residSS <- (deviance(max.vol))
#Defines the R^2 values for both models
R2 <- 1-(residSS/totalSS)
R2

#Max volume stomach capacity by length 
V <- 1.464e-10*length^4.204
wae <- cbind(wae, V)

#weight of fish+stomach contents (max weight) if the fish had filled its stomach to max volume
m.weight <- (weight-st_weight)+(V/1.05)
wae <- cbind(wae, m.weight)
rel_weight.max <- (m.weight/(10^(-5.453+3.180*(log10(length)))))*100
wae <- cbind(wae, rel_weight.max)

plotmeans(rel_weight.max~cls.5, na.rm=T, p=0.95, n.label=F, connect=F, gap=0, barcol="black")

################################################################################
#For SMB data
attach(stomach)
detach(wae)
smb <- stomach[species=="SMB",]
attach(smb)
detach(stomach)

#p.content <- 100*(st_weight/weight)
#smb <- cbind(smb, p.content)

rel_weight <- (weight/(10^(-5.329+3.2*(log10(length)))))*100
smb <- cbind(smb, rel_weight)

rel_weight_empty <- ((weight-st_weight)/(10^(-5.329+3.2*(log10(length)))))*100
smb <- cbind(smb, rel_weight_empty)

t.test(rel_weight, rel_weight_empty, paired=T)

psd <- with(smb,
            ifelse((length>=180)&(length<280), "S-Q",
                   ifelse((length>=280)&(length<350), "Q-P",
                          ifelse((length>=350)&(length<430), "P-M",
                                 ifelse((length>=430)&(length<510), "M-T",
                                        ifelse(length>=510, ">T",
                                               "substock"))))))
smb <- cbind(psd, smb)

tapply(rel_weight, psd, mean)
tapply(rel_weight_empty, psd, mean)

((tapply(rel_weight_empty, psd, mean)-tapply(rel_weight, psd, mean))/tapply(rel_weight, psd, mean))*100

#Function to determine sample size
values=function(x)
{
  return(sum(!is.na(x)))
}
tapply(rel_weight, psd, values)

shapiro.test(rel_weight)

with(wae[psd=="substock",], 
     kruskal.test(rel_weight, g=rel_weight_empty))
with(wae[psd=="S-Q",],
     kruskal.test(rel_weight, g=rel_weight_empty))
with(wae[psd=="Q-P",],
     kruskal.test(rel_weight, g=rel_weight_empty))
with(wae[psd=="P-M",],
     kruskal.test(rel_weight, g=rel_weight_empty))
with(wae[psd=="M-T",],
     kruskal.test(rel_weight, g=rel_weight_empty))

#length class code
cls.50 <- with(smb,
               ifelse((length>=50)&(length<100),"075",
                      ifelse((length>=100)&(length<150),"125",
                             ifelse((length>=150)&(length<200),"175",
                                    ifelse((length>=200)&(length<250),"225",
                                           ifelse((length>=250)&(length<300),"275",
                                                  ifelse((length>=300)&(length<350),"325",
                                                         ifelse((length>=350)&(length<400),"375",
                                                                ifelse((length>=400)&(length<450),"425",
                                                                       ifelse((length>=450)&(length<500),"475",
                                                                              ifelse((length>=500)&(length<550),"525",
                                                                                     ifelse((length>=550)&(length<600),"575",
                                                                                            ifelse((length>=600)&(length<650),"625",
                                                                                                   ifelse((length>=650)&(length<700),"675",
                                                                                                          ifelse((length>=700)&(length<750),"725",
                                                                                                                 "substock")))))))))))))))

smb <- cbind(smb, cls.50)

###Max volume based on 10 mm length classes
#cls.1 <- as.numeric(as.character(cut(length,seq(70,500,by=9.99),labels=seq(75,495,by=10),include.lowest=TRUE)))
#smb <- cbind(smb, cls.1)
#x.1 <- c(75,85,95,105,115,125,135,145,155,165,175,185,195,205,215,225,235,245,255,
#  265,275,285,295,305,315,325,335,345,355,365,375,385,395,405,415,425,435,445,455,465,475,485)
#y.1 <- c(0.1260,0.2100,.9765,0.4725,0.7665,0.7875,1.1865,1.6590,3.2655,4.1790,4.4100,
#  3.7695,6.3945,5.3970,5.1975,11.4660,7.8960,10.6890,6.9615,11.7285,12.2010,13.7970,
#  14.5950,17.2095,11.8965,14.3745,25.6200,14.5005,30.0510,31.5525,22.7850,19.2885,22.8900,
#  30.5445,32.6235,31.5525,21.0000,28.9275,42.8925,67.0005,54.1485,67.2315)
#plot(y.1~x.1, pch=19, ylab="Maximum stomach contents (ml)", xlab="length class (10 mm)", 
#  ylim=c(0,70), xlim=c(0, 500),yaxs="i", xaxs="i", las=1, bty="n")
#mod.1 <- nls(y.1~a.1*x.1^b.1, start=list(a.1=0.00000001, b.1=3.35), trace=T)
#  mod.2=seq(min(length), max(length), 0.01)
#  lines(mod.2, predict(mod.1, list(x.1 = mod.2)))   
#  text(200, 50, expression(V==3.99e-07*(length)^3.032))
#  text(200, 47, expression(R^2==0.87))
#summary(mod.1)
#R2
#  1-(deviance(mod.1)/((length(x.1)-1)*var(y.1)))

vol <- 1.05*st_weight
smb <- cbind(smb, vol)

tapply(vol, cls.50, max)
x <- c(75,125,175,225,275,325,375,425,475)
y <- c(0.9765,1.6590,6.3945,11.4660,14.5950,25.62,31.5525,32.6235,67.2315)
plot(y~x, pch=19, ylab="Maximum stomach contents (ml)", xlab="length class (50 mm)", 
     ylim=c(0,70), xlim=c(0, 500),yaxs="i", xaxs="i", las=1, bty="n")
max.vol <- nls(y~a*x^b, start=list(a=0.00000001, b=3.35), trace=T)
mod <- seq(min(length), max(length), 0.01)
lines(mod, predict(max.vol, list(x = mod)))   
text(200, 50, expression(V==4.73e-06*(length)^2.65))
text(200, 47, expression(R^2==0.94))
summary(max.vol)

#R2 value for the model
#Defines the total sum of squares for both models (totalSS=(n-1)*var(y))
totalSS <- (length(x)-1)*var(y)
#residSS=deviance(nls model)
residSS <- (deviance(max.vol))
#Defines the R^2 values for both models
1-(residSS/totalSS)



V <- 4.727e-06*length^2.654
smb <- cbind(smb, V)

#weight of fish+stomach contents (max weight) if the fish had filled its stomach to max volume
m.weight <- (weight-st_weight)+(V/1.05)
smb <- cbind(smb, m.weight)
rel_weight.max <- (m.weight/(10^(-5.329+3.2*(log10(length)))))*100
smb <- cbind(smb, rel_weight.max)

plotmeans(rel_weight.max~cls.50, na.rm=T, p=0.95, n.label=F, connect=F, gap=0, barcol="black")

t.test(rel_weight, rel_weight.max)
t.test(rel_weight, rel_weight_empty)
((tapply(rel_weight.max, psd, mean)-tapply(rel_weight, psd, mean))/tapply(rel_weight, psd, mean))*100