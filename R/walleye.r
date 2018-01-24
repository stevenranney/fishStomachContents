  stomach <- read.table("C:/Users/Steven Harris Ranney/Desktop/Manuscripts/Wr and Stomach Contents/data/stomach.txt", header=T)
  attach(stomach)
  names(stomach)

  wae <- stomach[Species=="WAE",]
  attach(wae)
  detach(stomach)

  wae$Lake <- as.factor(wae$Lake)

  #Calculate Wr from Murphy et al. 1990
  Wr <- (Weight/(10^(-5.453+3.18*(log10(Length)))))*100
  wae <- cbind(wae, Wr)

  #Calculate Wr-stomach contents from Murphy et al. 1990
  WrE <- ((Weight-St.weight)/(10^(-5.453+3.18*(log10(Length)))))*100
  wae <- cbind(wae, WrE)

  #Assign PSD values
  psd <- with(wae,
              ifelse((Length>=250)&(Length<380), "S-Q",
              ifelse((Length>=380)&(Length<510), "Q-P",
              ifelse((Length>=510)&(Length<630), "P-M",
              ifelse((Length>=630)&(Length<760), "M-T",
              ifelse(Length>=760, "M-T",
              "Substock"))))))
  wae <- cbind(psd, wae)

  #Convert stomach contents weigh to volume
  vol <- 1.05*St.weight
  wae <- cbind(wae, vol)

  ### Max volume and length by length category and population
  waeo <- wae[order(wae$Lake, wae$psd, wae$vol), ]
  wae.vol.length <- aggregate(waeo[c("Weight", "Length", "vol")], waeo[c("Lake", "psd")], tail, 1)
  #wae.vol.length <- read.table("wae.vol.length.txt")

  #max.vol.q <- nlrq(vol~a*Length^b, data=wae.vol.length, start=list(a=0.0000005, b=3.1), tau=0.75, trace=TRUE,
  #                   control=nls.control(maxiter=10000))

  #nonLinear regression of estimated stomach volume as a function of length
  max.vol.nls <- nls(vol~a*Length^b, data=wae.vol.length, start=list(a=0.00000005, b=3.0), trace=TRUE,
                     control=nls.control(maxiter=10000))

  with(wae.vol.length,
       plot(vol~Length, pch=19, ylab="Maximum stomach contents (ml)", xlab="Total length (mm)",
            ylim=c(0,200), xlim=c(200,700),yaxs="i", xaxs="i", las=1, bty="n"))
       mod <- seq(min(Length), max(Length), 0.01)
    #  lines(mod, predict(max.vol.q, list(Length = mod)))     
      lines(mod, predict(max.vol.nls, list(Length = mod)), lty=2)
  #    abline(v=250, lty=2)
  #    abline(v=380, lty=2)
  #    abline(v=510, lty=2)
  #    abline(v=630, lty=2)
  #  legend(225, 190, legend=c(expression(paste(3^rd, " Quantile")), "Nonlinear Least Squares"), lty=c(1,2))
  #  legend(225, 190, legend=c(expression(italic(Q[3])), "NLS"), lty=c(1,2))
    text(300, 190, label="B", cex=2)
  
  R2(wae.vol.length$Length, wae.vol.length$vol, max.vol.nls)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
  
  #Create a vector of total weights MINUS the stomach contents weight
  wae$deltaWeight <- wae$Weight-wae$St.weight
  
  #Quantile regression

  require(quantreg)
  #quantMod <- nlrq(wae.vol.length$vol~a*wae.vol.length$Length^b, start=list(a=.01, b=3.0), tau=.5)
  #quantMod <- nlrq(St.weight~a*Length^b, data=wae, start=list(a=.1, b=3.0), tau=.95)
  quantMod95 <- nlrq(St.weight~a*deltaWeight^b, data=wae, start=list(a=.1, b=3.0), tau=.95)
  #quantMod5 <- nlrq(St.weight~a*Weight^b, data=wae, start=list(a=.1, b=3.0), tau=.5)
  #nlsMod <- nls(St.weight~a*Weight^b, data=wae, start=list(a=.1, b=3.0))


  plot(wae$St.weight~wae$deltaWeight, pch=16, ylab="Stomach contents weight (g)", xlab="Total weight (g)",
       ylim=c(0,200), xlim=c(0,5000),yaxs="i", xaxs="i", las=1, bty="n")
  mod <- seq(min(wae$deltaWeight), max(wae$deltaWeight), length.out = length(wae$St.weight))
#  lines(mod, predict(max.vol.q, list(Length = mod)))     
  lines(mod, predict(quantMod95, list(deltaWeight= mod)), lty=2, lwd=2)
#  lines(mod, predict(quantMod5, list(Weight = mod)), lty=2, lwd=2, col="blue")
#  lines(mod, predict(nlsMod, list(Weight = mod)), lty=2, lwd=2, col="green")

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

