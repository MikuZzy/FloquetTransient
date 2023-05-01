# This file uses the .Rdata file from simulated data with varying median difference, and generates 
# the distribution plot.

setwd("~/R/novogene_trimmed/07_Additional")
load('07mediandiff.Rdata')

meddiff <- -newpar6-par_value[4]-par_value[7]
Type16f <- FNs_withsel$Type16f/1000
Nonconcordant <- FNs_withsel$other/1000
Othertypes <- 1-Type16f-Nonconcordant

# 800*600, 07_mediandiff.png
plot(meddiff,Type16f,type='l',lwd=1.5,col='red',lty=2,
     xlab='Difference in Median Baseline Expression',ylab='Fraction',xlim=c(0,7),ylim=c(0,1))
points(meddiff,Type16f,col='red',pch=19)
lines(meddiff,Nonconcordant,lwd=1.5,col='black',lty=2)
points(meddiff,Nonconcordant,col='black',pch=19)
lines(meddiff,Othertypes,lwd=1.5,col='blue',lty=2)
points(meddiff,Othertypes,col='blue',pch=19)
legend(x=c(5,7),y=c(0.7,0.5),legend=c("Model 16F","Other Models","Non Concordant"),
       col=c('red','blue','black'),lty=2,lwd=1,pch=19)
