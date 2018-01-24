stomach <- read.table("C:\\Users\\sranney\\Desktop\\Manuscripts\\Wr and Stomach Contents\\data\\stomach.txt", header=T)
attach(stomach)
names(stomach)

smb <- stomach[Species=="SMB",]
attach(smb)
detach(stomach)


Wr <- (Weight/(10^(-5.329+3.2*(log10(Length)))))*100
smb <- cbind(smb, Wr)

WrE <- ((Weight-St.weight)/(10^(-5.329+3.2*(log10(Length)))))*100
smb <- cbind(smb, WrE)

psd <- with(smb,
  ifelse((Length>=180)&(Length<280), "1",
  ifelse((Length>=280)&(Length<350), "2",
  ifelse((Length>=350)&(Length<430), "3",
  ifelse((Length>=430)&(Length<510), "4",
  ifelse(Length>=510, ">T",
  "0"))))))
smb <- cbind(psd, smb)

vol <- 1.05*St.weight
smb <- cbind(smb, vol)

### Max volume and length by length category and population
smbo <- smb[order(smb$Lake, smb$psd, smb$vol), ]
smb.vol.length <- write.table(aggregate(smbo[c("Length", "vol")], smbo[c("Lake", "psd")], tail, 1), file="smb.vol.length.txt")
smb.vol.length <- read.table("smb.vol.length.txt")

#max.vol.q <- nlrq(vol~a*Length^b, data=smb.vol.length, start=list(a=0.00000005, b=3.0), tau=0.75, trace=TRUE)
max.vol.nls <- nls(vol~a*Length^b, data=smb.vol.length, start=list(a=0.00000005, b=3.0), trace=TRUE)

with(smb.vol.length,
  plot(vol~Length, pch=19, ylab="Maximum stomach contents (ml)", xlab="Total length (mm)",
  ylim=c(0,80), xlim=c(100, 500),yaxs="i", xaxs="i", las=1, bty="n"))
  mod <- seq(min(Length), max(Length), 0.01)
#  lines(mod, predict(max.vol.q, list(Length = mod)))
  lines(mod, predict(max.vol.nls, list(Length = mod)), lty=2)
#    abline(v=180, lty=2)
#    abline(v=280, lty=2)
#    abline(v=350, lty=2)
#    abline(v=430, lty=2)
#    text(140, 75, "substock")
#    text(230, 75, "S-Q")
#    text(315, 75, "Q-P")
#    text(390, 75, "P-M")
#    text(470, 75, "M-T")

R2(smb.vol.length$Length, smb.vol.length$vol, max.vol.nls)


#Another method to calculate Coeffecient of determination
#SST <- sum((smb.vol.length$vol-mean(predict(max.vol)))^2)
#SSE <- sum((smb.vol.length$vol-(predict(max.vol)))^2)
#1-(SSE/SST)


#Calculate maximum stomach weight of individual fish
max.st.weight <- (2.727e-06*Length^2.696)/1.05
smb <- cbind(smb, max.st.weight)

#Calculate Wr.max of individual fish (i.e., total weight minus stomach weight
#plus estimated maximum stomach weight)
Wr.max <- ((Weight-St.weight+max.st.weight)/(10^(-5.329+3.2*(log10(Length))))*100)
smb <- cbind(smb, Wr.max)

#Percent difference between WrS and Wr
((tapply(WrE, psd, median)-tapply(Wr, psd, median))/tapply(Wr, psd, median))*100

#Function to determine sample size
values=function(x)
  {
  return(sum(!is.na(x)))
  }
tapply(Wr, psd, values)

shapiro.test(Wr)

with(smb[psd=="0",],
  wilcox.test(Wr, y=WrE))
with(smb[psd=="1",],
  wilcox.test(Wr, y=WrE))
with(smb[psd=="2",],
  wilcox.test(Wr, y=WrE))
with(smb[psd=="3",],
  wilcox.test(Wr, y=WrE))
with(smb[psd=="4",],
  wilcox.test(Wr, y=WrE))

  #Comparing WrE to WrMAX
with(smb[psd=="0",],
  wilcox.test(WrE, y=Wr.max))
with(smb[psd=="1",],
  wilcox.test(WrE, y=Wr.max))
with(smb[psd=="2",],
  wilcox.test(WrE, y=Wr.max))
with(smb[psd=="3",],
  wilcox.test(WrE, y=Wr.max))
with(smb[psd=="4",],
  wilcox.test(WrE, y=Wr.max))

  #Percent difference between Wr and Wr.max
((tapply(WrE, psd, median)-tapply(Wr.max, psd, median))/tapply(Wr.max, psd, median))*100

  #Comparing WrS to WrMAX
with(smb[psd=="0",],
  wilcox.test(Wr, y=Wr.max))
with(smb[psd=="1",],
  wilcox.test(Wr, y=Wr.max))
with(smb[psd=="2",],
  wilcox.test(Wr, y=Wr.max))
with(smb[psd=="3",],
  wilcox.test(Wr, y=Wr.max))
with(smb[psd=="4",],
  wilcox.test(Wr, y=Wr.max))

    #Percent difference between Wr and Wr.max
((tapply(Wr, psd, median)-tapply(Wr.max, psd, median))/tapply(Wr.max, psd, median))*100

###############################################################################################
#Barplot of the whole shebang!
#Adjusted for MEDIANS and SE, not MEAN and 95% CI.

z.smb <- matrix(c(tapply(Wr, psd, median),
              tapply(WrE, psd, median),
              tapply(Wr.max, psd, median)),
              nrow=3, ncol=5, byrow=T, dimnames=list(c("Wr", "WrE", "WrMAX"), c("0","1","2","3","4")))

barplot(z.smb, beside=T, space=c(0,1), ylim=c(85,115), xlab="Length Category", axis.lty=1, ylab=expression(paste("Relative weight   ", (italic(W[r])))),
        legend=F, xpd=F, names.arg=c("Substock","S-Q","Q-P","P-M","M-T"))#, density=c(0,0,6), angle=c(0,0,45), col="black")
  abline(h=85)
#  legend(12,110, legend=c(expression(italic(W[r]), italic(W[r])[E], italic(W[r])[MAX])), fill=gray.colors(3), bty="n", cex=2)

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

#matrices of 25th and 75th percentiles
z.q.l.smb <- matrix(c(quantile(Wr[psd==0], 0.25, na.rm=T), quantile(Wr[psd==1], 0.25, na.rm=T), quantile(Wr[psd==2], 0.25, na.rm=T), quantile(Wr[psd==3], 0.25, na.rm=T), quantile(Wr[psd==4], 0.25, na.rm=T),
                  quantile(WrE[psd==0], 0.25, na.rm=T), quantile(WrE[psd==1], 0.25, na.rm=T), quantile(WrE[psd==2], 0.25, na.rm=T), quantile(WrE[psd==3], 0.25, na.rm=T), quantile(WrE[psd==4], 0.25, na.rm=T),
                  quantile(Wr.max[psd==0], 0.25, na.rm=T), quantile(Wr.max[psd==1], 0.25, na.rm=T), quantile(Wr.max[psd==2], 0.25, na.rm=T), quantile(Wr.max[psd==3], 0.25, na.rm=T), quantile(Wr.max[psd==4], 0.25, na.rm=T)),
                  nrow=3, ncol=5, byrow=T, dimnames=list(c("Wr", "WrE", "WrMAX"), c("0","1","2","3","4"))
                  )

z.q.u.smb <- matrix(c(quantile(Wr[psd==0], 0.75, na.rm=T), quantile(Wr[psd==1], 0.75, na.rm=T), quantile(Wr[psd==2], 0.75, na.rm=T), quantile(Wr[psd==3], 0.75, na.rm=T), quantile(Wr[psd==4], 0.75, na.rm=T),
                  quantile(WrE[psd==0], 0.75, na.rm=T), quantile(WrE[psd==1], 0.75, na.rm=T), quantile(WrE[psd==2], 0.75, na.rm=T), quantile(WrE[psd==3], 0.75, na.rm=T), quantile(WrE[psd==4], 0.75, na.rm=T),
                  quantile(Wr.max[psd==0], 0.75, na.rm=T), quantile(Wr.max[psd==1], 0.75, na.rm=T), quantile(Wr.max[psd==2], 0.75, na.rm=T), quantile(Wr.max[psd==3], 0.75, na.rm=T), quantile(Wr.max[psd==4], 0.75, na.rm=T)),
                  nrow=3, ncol=5, byrow=T, dimnames=list(c("Wr", "WrE", "WrMAX"), c("0","1","2","3","4"))
                  )

  errbar((c(1.5:3.5, 5.5:7.5, 9.5:11.5, 13.5:15.5, 17.5:19.5)),
          z.smb,
          z.q.u.smb,
          z.q.l.smb,
          add=T, pch=" ")

                  
#detach(smb)
#rm(list=ls(all=TRUE))

#Least-squares regression on log10-transformed data
#smb.vol.length <- read.table("smb.vol.length.txt")
#mod.1 <- lm(log10(smb.vol.length$vol)~log10(smb.vol.length$Length))
#plot(log10(smb.vol.length$vol)~log10(smb.vol.length$Length))
#  abline(mod.1)
#summary(mod.1)
