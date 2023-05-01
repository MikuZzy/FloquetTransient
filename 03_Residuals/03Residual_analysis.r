# This file uses the full model fit of experimental data (with low expression gene filtered out), and generates
# histograms and scatter plots of ratios of residual MAD, fitted oscillation amplitude and median expression level
# on Day4/5 between two temperatures.

setwd("~/R/novogene_trimmed")
source("src_funcs.r")
load('01novogene_job.Rdata')
outdata_filtered <- read.csv(file = '01Result_filtered_fullparafit.csv', header = TRUE)

setwd("~/R/novogene_trimmed/03_Residuals")

res_analysis <- as.data.frame(matrix(NA,nrow=nrow(outdata_filtered),ncol=8))
colnames(res_analysis) <- c("geneid","gene_name","18C_Mad","18C_OscAmp","18C_MedExp",
                            "25C_Mad","25C_OscAmp","25C_MedExp")
rownames(res_analysis) <- outdata_filtered$X
res_analysis$geneid <- outdata_filtered$geneid
res_analysis$gene_name <- outdata_filtered$gene_name
indDay45_25C <- which((time_comb<0) & (time_comb>-48))
indDay45_18C <- which((time_comb<120) & (time_comb>72))

for (i in seq(nrow(outdata_filtered))) {
  j <- outdata_filtered$X[i]
  value <- as.vector(c(unlist(allgenes_data_25C[j,]),unlist(allgenes_data_18C[j,])))
  res <- value - func(as.vector(unlist(outdata_filtered[i,4:12])),time_comb)
  
  res_analysis$`18C_Mad`[i] <- mean(abs(res[indDay45_18C]))
  res_analysis$`25C_Mad`[i] <- mean(abs(res[indDay45_25C]))
  res_analysis$`18C_OscAmp`[i] <- sqrt((outdata_filtered$cos_25C[i]+outdata_filtered$cos_18C[i])^2 +
                                       (outdata_filtered$sin_25C[i]+outdata_filtered$sin_18C[i])^2)
  res_analysis$`25C_OscAmp`[i] <- sqrt(outdata_filtered$cos_25C[i]^2+outdata_filtered$sin_25C[i]^2)
  res_analysis$`18C_MedExp`[i] <- outdata_filtered$int_25C[i] - outdata_filtered$decaycos[i] - 
                                  outdata_filtered$decayint[i] - outdata_filtered$cos_18C[i]
  res_analysis$`25C_MedExp`[i] <- outdata_filtered$int_25C[i]
  message(i,"/",nrow(outdata_filtered))
}

# apply cutoff on log ratios
log_MadRatio <- sapply(res_analysis$`18C_Mad`/res_analysis$`25C_Mad`,
                       function(x) max(min(log(x),3),-3))
log_AmpRatio <- sapply(res_analysis$`18C_OscAmp`/res_analysis$`25C_OscAmp`,
                       function(x) max(min(log(x),5),-5))
log_ExpRatio <- sapply(res_analysis$`18C_MedExp`/res_analysis$`25C_MedExp`,
                       function(x) max(min(log(max(x,0)),3),-3))

cor(x=log_MadRatio,y=log_AmpRatio,method="pearson") # 0.3448967
cor(x=log_MadRatio,y=log_ExpRatio,method="pearson") # 0.4497940
t.test(na.omit(log(res_analysis$`18C_Mad`/res_analysis$`25C_Mad`))) # mu=0.3293562, p-value < 2.2e-16
t.test(na.omit(log(res_analysis$`18C_OscAmp`/res_analysis$`25C_OscAmp`))) # mu=0.3616574, p-value < 2.2e-16
t.test(na.omit(log(res_analysis$`18C_MedExp`/res_analysis$`25C_MedExp`))) # mu=0.07214271, p-value < 2.2e-16

# 750*550, 03_MADratio.png
hist(log_MadRatio, breaks=seq(-3,3,0.25),main="Log Day4/5 Residual MAD Ratio of 18C/25C",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)

# 750*550, 03_OscAmpratio.png
hist(log_AmpRatio, breaks=seq(-5,5,0.25),main="Log Oscillation Amplitude Ratio of 18C/25C",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)

# 750*550, 03_MBEratio.png
hist(log_ExpRatio, breaks=seq(-3,3,0.25),main="Log Median Baseline Expression Ratio of 18C/25C",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)

# 750*550, 03_scatter01.png
plot(log_MadRatio,log_AmpRatio,type="p",xlim=c(-3,3),ylim=c(-5,5),
     xlab="Log Day4/5 Residual MAD Ratio of 18C/25C",
     ylab="Log Oscillation Amplitude Ratio of 18C/25C",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)

# 750*550, 03_scatter02.png
plot(log_MadRatio,log_ExpRatio,type="p",xlim=c(-3,3),ylim=c(-3,3),
     xlab="Log Day4/5 Residual MAD Ratio of 18C/25C",
     ylab="Log Median Baseline Expression Ratio of 18C/25C",
     cex.axis=1.5,cex.lab=1.5,cex.main=1.5)

save.image(file = '03Residual_analysis.Rdata')
