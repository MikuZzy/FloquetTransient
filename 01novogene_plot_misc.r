# This file generates plots of proteasomes described in network section.

setwd("~/R/novogene_trimmed")
library("ggplot2")
library("egg")
source("src_funcs.r")
load("~/R/novogene_trimmed/01novogene_job.Rdata")
outdata_concordant <- read.csv(file = '01Result_concordant_details_repro.csv', header = TRUE)
vnames <- c("Rpn6","Prosalpha4","Prosbeta3","Rpn3","Prosbeta7","Clk","Rpn9","Rpt3","Rpt5","Rpn8","Usp14","Rpn7")
ind <- sapply(vnames,function(x,y) which(y==x),y=outdata$gene_name)
counter <- 1
g <- vector('list',length=length(ind))

setwd("~/R/novogene_trimmed/05_GRN")
for (i in ind) {
  value_18C <- unlist(allgenes_data_18C[i,])
  value_25C <- unlist(allgenes_data_25C[i,])
  temp_data <- data.frame(
    'value' = c(value_25C,value_18C),
    'time' = time_comb
  )
  g[[counter]] <- get_plots_nlrq_type_nolabel(par = as.vector(na.omit(unlist(outdata3[i,3:11]))),
                                              data18C = value_18C,
                                              data25C = value_25C,
                                              time1 = time_comb,
                                              time2 = fitted_time,
                                              time3 = seq(0,120,0.5),
                                              model = outdata2$model[i],
                                              type = outdata2$type[i],
                                              lambda = outdata2$TRUE_lambda[i])
  counter <- counter + 1
}

g_total <- ggarrange(plots=g,nrow=length(ind)/2,labels=vnames,label.args=list(gp=grid::gpar(font=4,cex=1.8)))
ggsave("05_ExpPattern.pdf",plot=g_total,width=12,height=12)


