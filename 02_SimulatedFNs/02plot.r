# This file uses the output .csv files from simulation, and generates histograms of method performance.

setwd("~/R/novogene_trimmed/02_SimulatedFNs")
library("ggplot2")
library("reshape2")

PER_nosel <- read.csv("02Result_PER_nosel.csv",header=TRUE)
PER_withsel <- read.csv("02Result_PER_withsel.csv",header=TRUE)

PER_plotdata <- data.frame(
  "Without Model Selection" = PER_nosel$Type16/1000,
  "With Model Selection" = (PER_withsel$Type16f+PER_withsel$Type16p)/1000
)
melt_PERplotdata <- melt(PER_plotdata)

PER_concordant <- data.frame(
  "Not Concordant" = PER_withsel$other,
  "Other Models" = 1000-PER_withsel$other-PER_withsel$Type16f-PER_withsel$Type14f,
  "Model14F" = PER_withsel$Type14f,
  "Model16F" = PER_withsel$Type16f
)
PER_concordant <- melt(PER_concordant)
PER_concordant$b <- rep(seq(0.5,5,0.5),4)
colnames(PER_concordant) <- c("Type","value","b")

melt_PERplotdata$b <- rep(seq(0.5,5,0.5),2)
PER_plotdata$Type16.in.concordant.classifications <- (PER_withsel$Type16f+PER_withsel$Type16p)/(1000-PER_withsel$other)
PER_plotdata$b <- seq(0.5,5,0.5)
PER_plotdata$concordant.classifications <- (1000-PER_withsel$other)/1000
PER_plotdata$Type16f.in.concordant.classified.Type16 <- PER_withsel$Type16f/(PER_withsel$Type16f+PER_withsel$Type16p)
PER_palette <- c("#33a02c","#b2df8a","#1f78b4","#a6cee3")

PER_MAD <- data.frame(
  "MAD" = 1.5,
  "ylevel" = 0.025
)

# 1000*600, 02plot_PER01
PER_ggp01 <- ggplot(data=PER_plotdata)+
  geom_bar(data=melt_PERplotdata,aes(x=b,y=value,fill=variable),stat="identity",
           width=.4,position = "dodge")+
  geom_line(data=PER_plotdata,aes(x=b,y=Type16.in.concordant.classifications,
                                  col="Type16 in concordant classifications"),size=2)+
  geom_point(data=PER_plotdata,aes(x=b,y=Type16.in.concordant.classifications),color="black",size=4)+
  geom_point(data=PER_MAD,aes(x=MAD,y=ylevel),shape=25,fill="red",color="red", size=4)+
  labs(title= "",x="Noise MAD",y="Fraction of Type16")+
  scale_y_continuous(limits=c(0,1),sec.axis=sec_axis(~.*1,name="Percentage"))+
  scale_color_manual(name=NULL,
                     breaks=c("Type16 in concordant classifications"),
                     values=c("Type16 in concordant classifications"="black"))+
  scale_fill_manual(name=NULL,
                    values=c("Without.Model.Selection"="grey", "With.Model.Selection"="black"))+
  theme(legend.position="bottom",
        legend.text=element_text(size=10))
ggsave("02plot_PER01.pdf",width=8,height=4,plot=PER_ggp01)

# 800*600, 02plot_PER02
PER_ggp02 <- ggplot(data=PER_concordant)+
  geom_bar(data=PER_concordant,aes(x=b,y=value,fill=Type),stat="identity",
           width=.4,position="fill")+
  geom_point(data=PER_MAD,aes(x=MAD,y=ylevel),shape=25,fill="red",color="red", size=4)+
  labs(title= "",x="Noise MAD",y="Fraction of Each Type")+
  scale_fill_manual(values=PER_palette)+
  theme(legend.position="bottom",
        legend.text=element_text(size=10))
ggsave("02plot_PER02.pdf",width=8,height=4,plot=PER_ggp02)

############################################
############################################
############################################

HSF_nosel <- read.csv("02Result_HSF_nosel.csv",header=TRUE)
HSF_withsel <- read.csv("02Result_HSF_withsel.csv",header=TRUE)

HSF_plotdata <- data.frame(
  "Without Model Selection" = HSF_nosel$Type16/1000,
  "With Model Selection" = (HSF_withsel$Type16f+HSF_withsel$Type16p)/1000
)
melt_HSFplotdata <- melt(HSF_plotdata)

HSF_concordant <- data.frame(
  "Not Concordant" = HSF_withsel$other,
  "Other Models" = 1000-HSF_withsel$other-HSF_withsel$Type02f-HSF_withsel$Type08f
                  -HSF_withsel$Type14f-HSF_withsel$Type16f,
  "Model02F" = HSF_withsel$Type02f,
  "Model08F" = HSF_withsel$Type08f,
  "Model14F" = HSF_withsel$Type14f,
  "Model16F" = HSF_withsel$Type16f
)
HSF_concordant <- melt(HSF_concordant)
HSF_concordant$b <- rep(seq(0.5,5,0.5),6)
colnames(HSF_concordant) <- c("Type","value","b")

melt_HSFplotdata$b <- rep(seq(0.5,5,0.5),2)
HSF_plotdata$Type16.in.concordant.classifications <- (HSF_withsel$Type16f+HSF_withsel$Type16p)/(1000-HSF_withsel$other)
HSF_plotdata$b <- seq(0.5,5,0.5)
HSF_plotdata$concordant.classifications <- (1000-HSF_withsel$other)/1000
HSF_plotdata$Type16f.in.concordant.classified.Type16 <- HSF_withsel$Type16f/(HSF_withsel$Type16f+HSF_withsel$Type16p)
HSF_palette <- c("#8c510a","#d8b365","#f6e8c3","#c7eae5","#5ab4ac","#01665e")

HSF_MAD <- data.frame(
  "MAD" = 2.5,
  "ylevel" = 0.025
)

# 1000*600, 02plot_HSF01
HSF_ggp01 <- ggplot(data=HSF_plotdata)+
  geom_bar(data=melt_HSFplotdata,aes(x=b,y=value,fill=variable),stat="identity",
           width=.4,position = "dodge")+
  geom_line(data=HSF_plotdata,aes(x=b,y=Type16.in.concordant.classifications,
                                  col="Type16 in concordant classifications"),size=2)+
  geom_point(data=HSF_plotdata,aes(x=b,y=Type16.in.concordant.classifications),color="black",size=4)+
  geom_point(data=HSF_MAD,aes(x=MAD,y=ylevel),shape=25,fill="red",color="red", size=4)+
  labs(title= "",x="Noise MAD",y="Fraction of Type16")+
  scale_y_continuous(limits=c(0,1),sec.axis=sec_axis(~.*1,name="Percentage"))+
  scale_color_manual(name=NULL,
                     breaks=c("Type16 in concordant classifications"),
                     values=c("Type16 in concordant classifications"="black"))+
  scale_fill_manual(name=NULL,
                    values=c("Without.Model.Selection"="grey", "With.Model.Selection"="black"))+
  theme(legend.position="bottom",
        legend.text=element_text(size=10))
ggsave("02plot_HSF01.pdf",width=8,height=4,plot=HSF_ggp01)

# 800*600, 02plot_HSF02
HSF_ggp02 <- ggplot(data=HSF_concordant)+
  geom_bar(data=HSF_concordant,aes(x=b,y=value,fill=Type),stat="identity",
           width=.4,position="fill")+
  geom_point(data=HSF_MAD,aes(x=MAD,y=ylevel),shape=25,fill="red",color="red", size=4)+
  labs(title= "",x="Noise MAD",y="Fraction of Each Type")+
  scale_fill_manual(values=HSF_palette)+
  theme(legend.position="bottom",
        legend.text=element_text(size=10))
ggsave("02plot_HSF02.pdf",width=8,height=4,plot=HSF_ggp02)

