# This file loads the .Rdata files run on HPC cluster, and output figures of distribution of lambdas 
# under each of the three samplin schemes.

setwd("~/R/novogene_trimmed/00_SyntheticData")
library(ggplot2)
library(egg)
get_counts <- function(x,breaks,omitNA=TRUE,cutoff=TRUE) {
  data <- sapply(x,min,breaks[length(breaks)])
  result <- rep(NA,length(breaks)-1)
  for (i in seq(length(breaks)-1)) {
    result[i] <- sum(data>breaks[i] & data<=breaks[i+1], na.rm=TRUE)
  }
  return(result)
}

br <- seq(0,1,0.05)
plot_br <- seq(0.025,0.995,0.05)

#####################################
#####################################
# plot of case no lambda
load("~/R/novogene_trimmed/00_SyntheticData/00.02SyntheticData_nolambda.Rdata")
df_CLK <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,1],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,1],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,1],breaks=br)
)
str_CLK <- c(paste('1hr,',sum(df_CLK$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_CLK$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_CLK$Exp),'/100',sep=''))
df_PER <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,2],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,2],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,2],breaks=br)
)
str_PER <- c(paste('1hr,',sum(df_PER$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_PER$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_PER$Exp),'/100',sep=''))
df_VRI <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,3],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,3],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,3],breaks=br)
)
str_VRI <- c(paste('1hr,',sum(df_VRI$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_VRI$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_VRI$Exp),'/100',sep=''))
df_PDP <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,4],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,4],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,4],breaks=br)
)
str_PDP <- c(paste('1hr,',sum(df_PDP$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_PDP$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_PDP$Exp),'/100',sep=''))

g1 <- ggplot(data=df_CLK,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_CLK[1]),size=3) + geom_line(aes(y=X1hr,colour=str_CLK[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_CLK[2]),size=3) + geom_line(aes(y=X2hr,colour=str_CLK[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_CLK[3]),size=3) + geom_line(aes(y=Exp,colour=str_CLK[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('CLK') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_CLK,values=c('black','blue','red'))

g2 <- ggplot(data=df_PER,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_PER[1]),size=3) + geom_line(aes(y=X1hr,colour=str_PER[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_PER[2]),size=3) + geom_line(aes(y=X2hr,colour=str_PER[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_PER[3]),size=3) + geom_line(aes(y=Exp,colour=str_PER[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('PER') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_PER,values=c('black','blue','red'))

g3 <- ggplot(data=df_VRI,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_VRI[1]),size=3) + geom_line(aes(y=X1hr,colour=str_VRI[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_VRI[2]),size=3) + geom_line(aes(y=X2hr,colour=str_VRI[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_VRI[3]),size=3) + geom_line(aes(y=Exp,colour=str_VRI[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('VRI') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_VRI,values=c('black','blue','red'))

g4 <- ggplot(data=df_PDP,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_PDP[1]),size=3) + geom_line(aes(y=X1hr,colour=str_PDP[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_PDP[2]),size=3) + geom_line(aes(y=X2hr,colour=str_PDP[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_PDP[3]),size=3) + geom_line(aes(y=Exp,colour=str_PDP[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('PDP') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_PDP,values=c('black','blue','red'))

g <- ggarrange(g1,g2,g3,g4,nrow=2)
ggsave(file='00.05Plot_nolambda.pdf',plot=g,height=10,width=16)

#####################################
#####################################
# plot of case lambda=0.50
load("~/R/novogene_trimmed/00_SyntheticData/00.03SyntheticData_lambda0.50.Rdata")
df_CLK <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,1],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,1],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,1],breaks=br)
)
str_CLK <- c(paste('1hr,',sum(df_CLK$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_CLK$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_CLK$Exp),'/100',sep=''))
df_PER <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,2],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,2],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,2],breaks=br)
)
str_PER <- c(paste('1hr,',sum(df_PER$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_PER$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_PER$Exp),'/100',sep=''))
df_VRI <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,3],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,3],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,3],breaks=br)
)
str_VRI <- c(paste('1hr,',sum(df_VRI$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_VRI$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_VRI$Exp),'/100',sep=''))
df_PDP <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,4],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,4],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,4],breaks=br)
)
str_PDP <- c(paste('1hr,',sum(df_PDP$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_PDP$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_PDP$Exp),'/100',sep=''))

g1 <- ggplot(data=df_CLK,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_CLK[1]),size=3) + geom_line(aes(y=X1hr,colour=str_CLK[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_CLK[2]),size=3) + geom_line(aes(y=X2hr,colour=str_CLK[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_CLK[3]),size=3) + geom_line(aes(y=Exp,colour=str_CLK[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('CLK') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_CLK,values=c('black','blue','red'))

g2 <- ggplot(data=df_PER,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_PER[1]),size=3) + geom_line(aes(y=X1hr,colour=str_PER[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_PER[2]),size=3) + geom_line(aes(y=X2hr,colour=str_PER[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_PER[3]),size=3) + geom_line(aes(y=Exp,colour=str_PER[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('PER') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_PER,values=c('black','blue','red'))

g3 <- ggplot(data=df_VRI,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_VRI[1]),size=3) + geom_line(aes(y=X1hr,colour=str_VRI[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_VRI[2]),size=3) + geom_line(aes(y=X2hr,colour=str_VRI[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_VRI[3]),size=3) + geom_line(aes(y=Exp,colour=str_VRI[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('VRI') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_VRI,values=c('black','blue','red'))

g4 <- ggplot(data=df_PDP,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_PDP[1]),size=3) + geom_line(aes(y=X1hr,colour=str_PDP[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_PDP[2]),size=3) + geom_line(aes(y=X2hr,colour=str_PDP[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_PDP[3]),size=3) + geom_line(aes(y=Exp,colour=str_PDP[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('PDP') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_PDP,values=c('black','blue','red'))

g <- ggarrange(g1,g2,g3,g4,nrow=2)
ggsave(file='00.05Plot_lambda0.50.pdf',plot=g,height=10,width=16)

#####################################
#####################################
# plot of case lambda=0.05
load("~/R/novogene_trimmed/00_SyntheticData/00.04SyntheticData_lambda0.05.Rdata")
df_CLK <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,1],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,1],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,1],breaks=br)
)
str_CLK <- c(paste('1hr,',sum(df_CLK$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_CLK$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_CLK$Exp),'/100',sep=''))
df_PER <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,2],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,2],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,2],breaks=br)
)
str_PER <- c(paste('1hr,',sum(df_PER$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_PER$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_PER$Exp),'/100',sep=''))
df_VRI <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,3],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,3],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,3],breaks=br)
)
str_VRI <- c(paste('1hr,',sum(df_VRI$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_VRI$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_VRI$Exp),'/100',sep=''))
df_PDP <- data.frame(
  'lambda' = plot_br,
  'X1hr' = get_counts(lambdasave_1hr[,4],breaks=br),
  'X2hr' = get_counts(lambdasave_2hr[,4],breaks=br),
  'Exp' = get_counts(lambdasave_exp[,4],breaks=br)
)
str_PDP <- c(paste('1hr,',sum(df_PDP$X1hr),'/100',sep=''),
             paste('2hr,',sum(df_PDP$X2hr),'/100',sep=''),
             paste('Exp,',sum(df_PDP$Exp),'/100',sep=''))

g1 <- ggplot(data=df_CLK,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_CLK[1]),size=3) + geom_line(aes(y=X1hr,colour=str_CLK[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_CLK[2]),size=3) + geom_line(aes(y=X2hr,colour=str_CLK[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_CLK[3]),size=3) + geom_line(aes(y=Exp,colour=str_CLK[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('CLK') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_CLK,values=c('black','blue','red'))

g2 <- ggplot(data=df_PER,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_PER[1]),size=3) + geom_line(aes(y=X1hr,colour=str_PER[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_PER[2]),size=3) + geom_line(aes(y=X2hr,colour=str_PER[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_PER[3]),size=3) + geom_line(aes(y=Exp,colour=str_PER[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('PER') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_PER,values=c('black','blue','red'))

g3 <- ggplot(data=df_VRI,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_VRI[1]),size=3) + geom_line(aes(y=X1hr,colour=str_VRI[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_VRI[2]),size=3) + geom_line(aes(y=X2hr,colour=str_VRI[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_VRI[3]),size=3) + geom_line(aes(y=Exp,colour=str_VRI[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('VRI') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_VRI,values=c('black','blue','red'))

g4 <- ggplot(data=df_PDP,aes(x=lambda)) + 
  geom_point(aes(y=X1hr,colour=str_PDP[1]),size=3) + geom_line(aes(y=X1hr,colour=str_PDP[1]),lwd=1.5) + 
  geom_point(aes(y=X2hr,colour=str_PDP[2]),size=3) + geom_line(aes(y=X2hr,colour=str_PDP[2]),lwd=1.5) + 
  geom_point(aes(y=Exp,colour=str_PDP[3]),size=3) + geom_line(aes(y=Exp,colour=str_PDP[3]),lwd=1.5) +
  scale_y_continuous(breaks=seq(0,100,20)) + 
  xlab('lambda') + ylab('counts') + ggtitle('PDP') +
  coord_cartesian(xlim=c(0,1),ylim=c(0,100)) +
  scale_colour_manual("",breaks=str_PDP,values=c('black','blue','red'))

g <- ggarrange(g1,g2,g3,g4,nrow=2)
ggsave(file='00.05Plot_lambda0.05.pdf',plot=g,height=10,width=16)










