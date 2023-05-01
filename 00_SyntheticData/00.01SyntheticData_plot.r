# This file plots the synthetic data from Matlab with three different changes: transient, fast and slow.

setwd("~/R/novogene_trimmed")
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.0",.libPaths()))
library("ggplot2")
library("egg")

setwd("~/R/novogene_trimmed/00_SyntheticData")

SynData <- read.csv(file='Data_Synthetic_nolambda.csv',header=TRUE)
SynData[,2:5] <- log(SynData[,2:5]+1)
g1 <- ggplot(data=SynData,aes(time,CLK)) + 
  geom_point(size=2,col="black") + geom_line(lwd=1.5,col="black") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g2 <- ggplot(data=SynData,aes(time,PER_total)) + 
  geom_point(size=2,col="blue") + geom_line(lwd=1.5,col="blue") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g3 <- ggplot(data=SynData,aes(time,VRI)) + 
  geom_point(size=2,col="red") + geom_line(lwd=1.5,col="red") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g4 <- ggplot(data=SynData,aes(time,PDP)) + 
  geom_point(size=2,col="cyan") + geom_line(lwd=1.5,col="cyan") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g <- ggarrange(plots=list(g1,g2,g3,g4),nrow=2,ncol=2)
ggsave(filename="SynData_nolambda_raw.pdf",plot=g,width=12,height=6)


SynData <- read.csv(file='Data_Synthetic_lambda0.50.csv',header=TRUE)
SynData[,2:5] <- log(SynData[,2:5]+1)
g1 <- ggplot(data=SynData,aes(time,CLK)) + 
  geom_point(size=2,col="black") + geom_line(lwd=1.5,col="black") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g2 <- ggplot(data=SynData,aes(time,PER_total)) + 
  geom_point(size=2,col="blue") + geom_line(lwd=1.5,col="blue") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g3 <- ggplot(data=SynData,aes(time,VRI)) + 
  geom_point(size=2,col="red") + geom_line(lwd=1.5,col="red") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g4 <- ggplot(data=SynData,aes(time,PDP)) + 
  geom_point(size=2,col="cyan") + geom_line(lwd=1.5,col="cyan") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g <- ggarrange(plots=list(g1,g2,g3,g4),nrow=2,ncol=2)
ggsave(filename="SynData_lambda0.50_raw.pdf",plot=g,width=12,height=6)


SynData <- read.csv(file='Data_Synthetic_lambda0.05.csv',header=TRUE)
SynData[,2:5] <- log(SynData[,2:5]+1)
g1 <- ggplot(data=SynData,aes(time,CLK)) + 
  geom_point(size=2,col="black") + geom_line(lwd=1.5,col="black") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g2 <- ggplot(data=SynData,aes(time,PER_total)) + 
  geom_point(size=2,col="blue") + geom_line(lwd=1.5,col="blue") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g3 <- ggplot(data=SynData,aes(time,VRI)) + 
  geom_point(size=2,col="red") + geom_line(lwd=1.5,col="red") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g4 <- ggplot(data=SynData,aes(time,PDP)) + 
  geom_point(size=2,col="cyan") + geom_line(lwd=1.5,col="cyan") +
  scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
        panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g <- ggarrange(plots=list(g1,g2,g3,g4),nrow=2,ncol=2)
ggsave(filename="SynData_lambda0.05_raw.pdf",plot=g,width=12,height=6)








