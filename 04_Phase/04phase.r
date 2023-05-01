# This file generates heat maps and phase scatter plot of concordant cycling genes.

library('ggplot2')
setwd("~/R/novogene_trimmed/00_Reproducibility")
result_filter <- read.csv(file="00ReproducibleGenes.csv")

setwd("~/R/novogene_trimmed")
source("src_funcs.r")
results <- read.csv(file = '01Result_filtered_fullparafit.csv', header = TRUE)
results2 <- read.csv(file = '01Result_filtered_modelAICc.csv', header = TRUE)
novogene2021 <- read.table(file = "Data_novogene2021_trimmed_dedup_tpm.tsv", header = TRUE)
gene_keys <- read.csv(file = "Data_SampleKey_v3.csv", header = TRUE)

# 18C Data from current data set
sample_index_18C <- which(gene_keys$Temperature == '18C')
time_18C <- (gene_keys$Day[sample_index_18C]-1)*24 + gene_keys$Hour[sample_index_18C]
# 25C Data from current data set
sample_index_25C <- which(gene_keys$Temperature == '25C')
time_25C <- (gene_keys$Day[sample_index_25C]-1)*24 + gene_keys$Hour[sample_index_25C]

# ZT for all data
time_comb <- sort(unique(c(time_25C-120,time_18C)))
n <- length(time_comb)

setwd("~/R/novogene_trimmed/04_Phase")

# get concordant cycling genes
sig <- 0.05
ind <- which(results$pval_cycle_25C <= sig &
             results$pval_cycle_18C <= sig &
             sapply(results2$model,function(x) as.numeric(substr(x,1,2))) == results2$type &
             apply(results2[,5:24],1,function(x) diff(sort(x))[1]) > 2 &
             results$geneid %in% result_filter$geneid)
DC <- results$pval_cycle_diff[ind] <= sig

temp_results <- results[ind,]
temp_data <- novogene2021[temp_results$X,]

ngene <- nrow(temp_results)
PRCdata <- matrix(NA, nrow=ngene, ncol=4)
colnames(PRCdata) <- c('geneid','Phase','PhaseRight','Response')
PRCdata <- as.data.frame(PRCdata)
PRCdata$geneid <- temp_data$geneid

# fitted phase
for (i in seq(ngene)) {
  PRCdata[i,2:4] <- get_PRC(temp_results$cos_25C[i],temp_results$sin_25C[i],
                         temp_results$cos_18C[i],temp_results$sin_18C[i])
}
PRCdata$Phase_proj <- PRCdata$Phase + 24
PRCdata$PhaseRight_proj <- PRCdata$PhaseRight + 24
PRCdata$color <- ifelse(DC,"D","ND")
PRCdata <- PRCdata[order(PRCdata$PhaseRight),]
PRCdata$phaseorder <- seq(nrow(PRCdata),1)
PRCdata <- PRCdata[order(PRCdata$geneid),]

coreclock <- which(PRCdata$geneid %in% c('FBgn0023076','FBgn0016076','FBgn0003068','FBgn0014396'))
PRCdata[coreclock,]
PRCdata$Phase[coreclock] - PRCdata$PhaseRight[coreclock]

heatmap_data <- as.data.frame(matrix(NA,nrow=length(ind)*n,ncol=4))
colnames(heatmap_data) <- c("geneid","ZT","TPM","Temp")
for (i in seq(length(ind))) {
  for (j in seq(n)) {
    y <- paste("00",as.character(PRCdata$phaseorder[i]),sep="")
    y <- substr(y,nchar(y)-2,nchar(y))
    heatmap_data[(i-1)*n+j,1] <- paste(y,PRCdata$geneid[i],sep="_")
    if (time_comb[j]<0) {
      x <- paste("00",as.character(time_comb[j]+120),sep="")
      x <- paste("ZT",substr(x,nchar(x)-2,nchar(x)),sep="")
      ind_time <- which(gene_keys$Temperature=='25C' & gene_keys$ZT.1==x)
      a <- paste("00",as.character(time_comb[j]+120),sep="")
      heatmap_data[(i-1)*n+j,4] <- "25C"
    } else {
      x <- paste("00",as.character(time_comb[j]),sep="")
      x <- paste("ZT",substr(x,nchar(x)-2,nchar(x)),sep="")
      ind_time <- which(gene_keys$Temperature=='18C' & 
                          gene_keys$ZT.1 == x)
      a <- paste("00",as.character(time_comb[j]),sep="")
      heatmap_data[(i-1)*n+j,4] <- "18C"
    }
    heatmap_data[(i-1)*n+j,2] <- substr(a,nchar(a)-2,nchar(a))
    heatmap_data[(i-1)*n+j,3] <- mean(as.numeric(temp_data[i,ind_time+2]))
    message(i,'/',length(ind),', ',j,'/',n)
  }
}
heatmap_data_25C <- heatmap_data[which(heatmap_data$Temp=="25C"),]
heatmap_data_18Cday12 <- heatmap_data[which(heatmap_data$Temp=="18C" &
                                              as.numeric(heatmap_data$ZT)<=48),]
heatmap_data_18Cday45 <- heatmap_data[which(heatmap_data$Temp=="18C" &
                                              as.numeric(heatmap_data$ZT)>=72),]
geneids <- sort(unique(heatmap_data$geneid))
for (i in seq(length(ind))) {
  geneid <- geneids[i]
  heatmap_data_25C$TPM[which(heatmap_data_25C$geneid==geneid)] <- 
    (heatmap_data_25C$TPM[which(heatmap_data_25C$geneid==geneid)] - mean(heatmap_data_25C$TPM[which(heatmap_data_25C$geneid==geneid)]))/sd(heatmap_data_25C$TPM[which(heatmap_data_25C$geneid==geneid)])
  heatmap_data_18Cday12$TPM[which(heatmap_data_18Cday12$geneid==geneid)] <- 
    (heatmap_data_18Cday12$TPM[which(heatmap_data_18Cday12$geneid==geneid)] - mean(heatmap_data_18Cday12$TPM[which(heatmap_data_18Cday12$geneid==geneid)]))/sd(heatmap_data_18Cday12$TPM[which(heatmap_data_18Cday12$geneid==geneid)])
  heatmap_data_18Cday45$TPM[which(heatmap_data_18Cday45$geneid==geneid)] <- 
    (heatmap_data_18Cday45$TPM[which(heatmap_data_18Cday45$geneid==geneid)] - mean(heatmap_data_18Cday45$TPM[which(heatmap_data_18Cday45$geneid==geneid)]))/sd(heatmap_data_18Cday45$TPM[which(heatmap_data_18Cday45$geneid==geneid)])
}

colnames(heatmap_data_25C) <- c("geneid","ZT","TPM_zscore","Temp")
colnames(heatmap_data_18Cday12) <- c("geneid","ZT","TPM_zscore","Temp")
colnames(heatmap_data_18Cday45) <- c("geneid","ZT","TPM_zscore","Temp")

# 850*800, 04_heatmap_25C.png
g1 <- ggplot(heatmap_data_25C, aes(ZT, geneid,fill=TPM_zscore)) + geom_tile() +
  scale_fill_viridis_c() + ylab("gene") + ggtitle("25C") +
  theme(axis.text.y=element_blank(),plot.title = element_text(hjust = 0.5))
ggsave("04_heatmap_25C.pdf",g1,width=8.5,height=8)

# 850*800, 04_heatmap_18Cday12.png
g2 <- ggplot(heatmap_data_18Cday12, aes(ZT, geneid,fill=TPM_zscore)) + geom_tile() +
  scale_fill_viridis_c() + ylab("gene") + ggtitle("18C Day1/2") +
  theme(axis.text.y=element_blank(),plot.title = element_text(hjust = 0.5))
ggsave("04_heatmap_18Cday12.pdf",g2,width=8.5,height=8)

# 850*800, 04_heatmap_18Cday45.png
g3 <- ggplot(heatmap_data_18Cday45, aes(ZT, geneid,fill=TPM_zscore)) + geom_tile() +
  scale_fill_viridis_c() + ylab("gene") + ggtitle("18C Day4/5") +
  theme(axis.text.y=element_blank(),plot.title = element_text(hjust = 0.5))
ggsave("04_heatmap_18Cday45.pdf",g3,width=8.5,height=8)

data_rug <- data.frame(
  "PhaseLeft" = c(PRCdata$Phase,PRCdata$Phase_proj,PRCdata$Phase,PRCdata$Phase_proj),
  "PhaseRight" = c(PRCdata$PhaseRight,PRCdata$PhaseRight,PRCdata$PhaseRight_proj,PRCdata$PhaseRight_proj),
  "color" = rep(PRCdata$color,4)
)

# 1000*800, 04_scatter.png
g4 <- ggplot(data=data_rug,aes(PhaseLeft,PhaseRight,color=color)) +
  geom_vline(xintercept=24,color="black",size=0.25) +
  geom_hline(yintercept=24,color="black",size=0.25) +
  geom_point() + xlab("Phase 25C") + ylab("Phase 18C") +
  scale_x_continuous(breaks=seq(0,48,6),minor_breaks=24) +
  scale_y_continuous(breaks=seq(0,48,6),minor_breaks=24) +
  theme(legend.position="bottom",legend.spacing.x = unit(0.5,'cm')) + 
  geom_rug(aes(color=color)) +
  scale_colour_manual(name="Concordant Genes",values=c("D"="black","ND"="red"),
                      labels=c(paste("Cycling and Differentially Cycling, ",sum(DC),"/",length(ind),sep=""),
                               paste("Cycling but not Differentially Cycling, ",sum(!DC),"/",length(ind),sep="")))
ggsave("04_scatter.pdf",g4,width=10,height=8)

colnames(PRCdata)[2] <- "Phase25C"
colnames(PRCdata)[3] <- "Phase18C"
write.csv(PRCdata,file='04_Result_phase.csv',row.names=FALSE)



