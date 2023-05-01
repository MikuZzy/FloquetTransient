# This file performs full model fitting on synthetic data curves with different residual MAD without doing
# any statistical inferences. Plots of fitted distribution of lambdas are generated.

setwd("~/R/novogene_trimmed")
library('quantreg')
library('parallel')
source('src_funcs.r')

gene_keys <- read.csv(file = 'Data_SampleKey_v3.csv', header = TRUE)
sample_index_18C <- which(gene_keys$Temperature == '18C')
sample_index_25C <- which(gene_keys$Temperature == '25C')
time_18C <- (gene_keys$Day[sample_index_18C]-1)*24 + gene_keys$Hour[sample_index_18C]
time_25C <- (gene_keys$Day[sample_index_25C]-1)*24 + gene_keys$Hour[sample_index_25C]
time_comb <- c(time_25C-120,time_18C)
n <- length(time_comb)

setwd("~/R/novogene_trimmed/06_MultiLambda")
lambdasave <- rep(NA,1000)
value_25C <- rep(3,length(time_25C))
# value_18C <- 1 + exp(-time_18C*0.5) + exp(-time_18C*0.05)
value_18C <- 1 + exp(-time_18C*0.4) + exp(-time_18C*0.04)
nretry <- 5
fluc <- 0.01
set.seed(0)

for (i in seq(1000)) {
  res <- rexp(n,rate=10)*sample(c(1,-1),n,replace=TRUE)
  temp_data <- data.frame(
    'value' = c(value_25C,value_18C) + res,
    'time' = time_comb
  )
  par_est <- get_guess(value = temp_data$value, time = temp_data$time,
                       coarse = seq(0.04,1,0.03), fine = seq(0.04,1,0.01))
  
  ind <- 0
  errmsg <- TRUE
  while (errmsg & ind < nretry) {
    ind <- ind + 1
    par_init <- par_est * runif(9, min = 1-fluc, max = 1+fluc)
    nlrqmodel <- optim(par = par_init, fn = obj, gr = grad, method = 'L-BFGS-B',
                       value = temp_data$value, time = temp_data$time,
                       lower = c(rep(NA,8),0.03), upper = c(rep(NA,8),1.5))
    if (nlrqmodel$convergence == 0) {
      errmsg <- FALSE
    }
  }
  lambdasave[i] <- nlrqmodel$par[9]
  message(i,"/1000")
}
max(lambdasave)
min(lambdasave)
# 750*550, 06_MAD0.1.png
hist(lambdasave,breaks=seq(0,0.3,0.01),main="MAD = 0.1",xlab="fitted lambda")



lambdasave2 <- rep(NA,1000)
for (i in seq(1000)) {
  res <- rexp(n,rate=2)*sample(c(1,-1),n,replace=TRUE)
  temp_data <- data.frame(
    'value' = c(value_25C,value_18C) + res,
    'time' = time_comb
  )
  par_est <- get_guess(value = temp_data$value, time = temp_data$time,
                       coarse = seq(0.04,1,0.03), fine = seq(0.04,1,0.01))
  
  ind <- 0
  errmsg <- TRUE
  while (errmsg & ind < nretry) {
    ind <- ind + 1
    par_init <- par_est * runif(9, min = 1-fluc, max = 1+fluc)
    nlrqmodel <- optim(par = par_init, fn = obj, gr = grad, method = 'L-BFGS-B',
                       value = temp_data$value, time = temp_data$time,
                       lower = c(rep(NA,8),0.03), upper = c(rep(NA,8),1.5))
    if (nlrqmodel$convergence == 0) {
      errmsg <- FALSE
    }
  }
  lambdasave2[i] <- nlrqmodel$par[9]
  message(i,"/1000")
}
max(lambdasave2)
min(lambdasave2)

# 750*550, 06_MAD0.5.png
hist(lambdasave2,breaks=seq(0,1.1,0.01),main="MAD = 0.5",xlab="fitted lambda")

save.image(file="06MultiLambda.Rdata")


