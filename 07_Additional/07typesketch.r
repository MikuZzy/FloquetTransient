# This file generates plots of each candidate model used in model selection, only for illustration purposes.

setwd("~/R/novogene_trimmed")
library("ggplot2")
library("egg")
source("src_funcs.r")

setwd("~/R/novogene_trimmed/07_Additional")
time_comb <- seq(-48,48,0.5)

g01 <- ggplot(data=data.frame("time"=time_comb,"value"=func_fullmodel_type(c(1,1,1,0.2),time_comb,1)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#332288',lwd=2) + coord_cartesian(ylim=c(0.5,1.5)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  # theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
  #       panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
  #       panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
g02 <- ggplot(data=data.frame("time"=time_comb,"value"=func_fullmodel_type(c(1,-2,1,1,0.2),time_comb,2)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#88CCEE',lwd=2) + coord_cartesian(ylim=c(0.5,3)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g07 <- ggplot(data=data.frame("time"=time_comb,"value"=func_fullmodel_type(c(1,-9,2,9,1,0.3),time_comb,7)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#44AA99',lwd=2) + coord_cartesian(ylim=c(-2.5,5)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g08 <- ggplot(data=data.frame("time"=time_comb,"value"=func_fullmodel_type(c(1,-2,-6,2,-1.5,1,0.2),time_comb,8)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#117733',lwd=2) + coord_cartesian(ylim=c(0,5)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g11 <- ggplot(data=data.frame("time"=time_comb,"value"=func_fullmodel_type(c(1,0,-1,1,1,0.3),time_comb,11)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#999933',lwd=2) + coord_cartesian(ylim=c(-2,4)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g12 <- ggplot(data=data.frame("time"=time_comb,"value"=func_fullmodel_type(c(1,-2,-1,-1,1,1,0.3),time_comb,12)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#DDCC77',lwd=2) + coord_cartesian(ylim=c(-2,4)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g13 <- ggplot(data=data.frame("time"=time_comb,"value"=func_fullmodel_type(c(1,-1,1,-6,2,0.3),time_comb,13)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#CC6677',lwd=2) + coord_cartesian(ylim=c(-2,4.5)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g14 <- ggplot(data=data.frame("time"=time_comb,"value"=func_fullmodel_type(c(1,-1,1,-4,2,2,0.3),time_comb,14)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#882255',lwd=2) + coord_cartesian(ylim=c(-1,6)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g15 <- ggplot(data=data.frame("time"=time_comb,"value"=func_fullmodel_type(c(1,-1,1,-4,2,2,1,0.5),time_comb,15)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#AA4499',lwd=2) + coord_cartesian(ylim=c(-2,5)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g16 <- ggplot(data=data.frame("time"=time_comb,"value"=func_fullmodel_type(c(1,-1,1,-4,2,4,1,1,0.3),time_comb,16)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#DDDDDD',lwd=2) + coord_cartesian(ylim=c(-3,4)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

g <- ggarrange(g01,g02,g07,g08,g11,g12,g13,g14,g15,g16,nrow=2,byrow=FALSE,
               labels=c('Type 01','Type 02','Type 07','Type 08','Type 11','Type 12','Type 13','Type 14','Type 15','Type 16'),
               label.args=list(gp=grid::gpar(font=2,cex=1,fontsize=8)))
ggsave('07_typesketch01_old_nogrid.pdf',plot=g,width=10,height=4)
# https://personal.sron.nl/~pault/#sec:qualitative

g01_p <- ggplot(data=data.frame("time"=time_comb,"value"=func_partialmodel_type(c(1),time_comb,1)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#332288',lwd=2) + coord_cartesian(ylim=c(0.5,1.5)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g02_p <- ggplot(data=data.frame("time"=time_comb,"value"=func_partialmodel_type(c(1,1),time_comb,2)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#88CCEE',lwd=2) + coord_cartesian(ylim=c(0.5,3)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g07_p <- ggplot(data=data.frame("time"=time_comb,"value"=func_partialmodel_type(c(1,1),time_comb,7)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#44AA99',lwd=2) + coord_cartesian(ylim=c(-2.5,5)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g08_p <- ggplot(data=data.frame("time"=time_comb,"value"=func_partialmodel_type(c(1,-1.5,1),time_comb,8)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#117733',lwd=2) + coord_cartesian(ylim=c(0,5)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g11_p <- ggplot(data=data.frame("time"=time_comb,"value"=func_partialmodel_type(c(1,-1),time_comb,11)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#999933',lwd=2) + coord_cartesian(ylim=c(-2,4)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g12_p <- ggplot(data=data.frame("time"=time_comb,"value"=func_partialmodel_type(c(1,-2,-1),time_comb,12)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#DDCC77',lwd=2) + coord_cartesian(ylim=c(-2,4)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL,y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g13_p <- ggplot(data=data.frame("time"=time_comb,"value"=func_partialmodel_type(c(1,-1,1),time_comb,13)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#CC6677',lwd=2) + coord_cartesian(ylim=c(-2,4.5)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g14_p <- ggplot(data=data.frame("time"=time_comb,"value"=func_partialmodel_type(c(1,-1,1,2),time_comb,14)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#882255',lwd=2) + coord_cartesian(ylim=c(-1,6)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g15_p <- ggplot(data=data.frame("time"=time_comb,"value"=func_partialmodel_type(c(1,-1,1,1),time_comb,15)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#AA4499',lwd=2) + coord_cartesian(ylim=c(-2,5)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g16_p <- ggplot(data=data.frame("time"=time_comb,"value"=func_partialmodel_type(c(1,-1,1,1,1),time_comb,16)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#DDDDDD',lwd=2) + coord_cartesian(ylim=c(-3,4)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(x=NULL, y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

g <- ggarrange(g01_p,g01,g02_p,g02,g07_p,g07,g08_p,g08,g11_p,g11,
               g12_p,g12,g13_p,g13,g14_p,g14,g15_p,g15,g16_p,g16,nrow=4,byrow=FALSE,
               labels=c('Model 01P','Model 01F','Model 02P','Model 02F','Model 07P','Model 07F',
                        'Model 08P','Model 08F','Model 11P','Model 11F','Model 12P','Model 12F',
                        'Model 13P','Model 13F','Model 14P','Model 14F','Model 15P','Model 15F',
                        'Model 16P','Model 16F'),
               label.args=list(gp=grid::gpar(font=2,cex=1,fontsize=8)))
ggsave('07_typesketch02_old_nogrid.pdf',plot=g,width=10,height=8)

g07 <- ggplot(data=data.frame("time"=time_comb,"value"=func_fullmodel_type(c(1,-9,2,9,1,0.3),time_comb,7)),aes(time,value)) +
  geom_vline(xintercept=0,colour="grey") + 
  geom_line(color='#44AA99',lwd=2) + coord_cartesian(ylim=c(-2.5,5)) + scale_x_continuous(breaks=seq(-48,48,by=24)) +
  labs(y=NULL) + 
  theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

g <- ggarrange(g01_p,g01,g07_p,g07,g13_p,g13,nrow=2,byrow=FALSE,
               labels=c('Model 01P','Model 01F','Model 07P','Model 07F','Model 13P','Model 13F'),
               label.args=list(gp=grid::gpar(font=2,cex=1,fontsize=8)))
ggsave('07_typesketch03_old_nogrid.pdf',plot=g,width=6,height=4)



