# This file uses the output from main script and generates plots of individual genes.

setwd("~/R/novogene_trimmed")
library("ggplot2")
library("egg")
source("src_funcs.r")
load("~/R/novogene_trimmed/01novogene_job.Rdata")
outdata_concordant <- read.csv(file = '01Result_concordant_details_repro.csv', header = TRUE)
ind <- which(outdata$gene_name %in% c('Clk','vri','per','tim'))

setwd("~/R/novogene_trimmed/01_Plots")

for (i in seq(ngene)) {
  
  if (!outdata$small_tpm[i]) {
    value_18C <- unlist(allgenes_data_18C[i,])
    value_25C <- unlist(allgenes_data_25C[i,])
    temp_data <- data.frame(
      'value' = c(value_25C,value_18C),
      'time' = time_comb
    )
    
    # plots
    g1 <- get_plots_nlrq(result = outdata[i,],
                         data18C = value_18C,
                         data25C = value_25C,
                         time1 = time_comb,
                         time2 = fitted_time,
                         time3 = seq(0,120,0.5))
    g2 <- get_plots_nlrq_type(par = as.vector(na.omit(unlist(outdata3[i,3:11]))),
                              data18C = value_18C,
                              data25C = value_25C,
                              time1 = time_comb,
                              time2 = fitted_time,
                              time3 = seq(0,120,0.5),
                              model = outdata2$model[i],
                              type = outdata2$type[i],
                              lambda = outdata2$TRUE_lambda[i])
    g <- ggarrange(g1,g2,nrow=2)
    filename <- paste('00000',i,'.pdf',sep='')
    filename <- substr(filename,start=nchar(filename)-8,stop=nchar(filename))
    ggsave(filename,plot=g,width=6,height=6)
    
    if (i %in% outdata_concordant$X) {
      setwd("~/R/novogene_trimmed/01_Plots_Concordant")
      g1 <- get_plots_nlrq_nolabel(result = outdata[i,],
                                   data18C = value_18C,
                                   data25C = value_25C,
                                   time1 = time_comb,
                                   time2 = fitted_time,
                                   time3 = seq(0,120,0.5))
      g2 <- get_plots_nlrq_type_nolabel(par = as.vector(na.omit(unlist(outdata3[i,3:11]))),
                                        data18C = value_18C,
                                        data25C = value_25C,
                                        time1 = time_comb,
                                        time2 = fitted_time,
                                        time3 = seq(0,120,0.5),
                                        model = outdata2$model[i],
                                        type = outdata2$type[i],
                                        lambda = outdata2$TRUE_lambda[i])
      g <- ggarrange(g1,g2,nrow=2)
      ggsave(filename,plot=g,width=6,height=6)
      setwd("~/R/novogene_trimmed/01_Plots")
    }
    
    if (i %in% ind) {
      setwd("~/R/novogene_trimmed/01_Plots_Coreclock")
      g1 <- get_plots_nlrq_nolabel(result = outdata[i,],
                                   data18C = value_18C,
                                   data25C = value_25C,
                                   time1 = time_comb,
                                   time2 = fitted_time,
                                   time3 = seq(0,120,0.5))
      ggsave(filename,plot=g1,width=6,height=3)
      setwd("~/R/novogene_trimmed/01_Plots")
    }
    
  }
  
  message(i)
  
}
