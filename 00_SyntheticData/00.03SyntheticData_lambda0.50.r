# This file works with synthetic data with fast changes, and compare three different sampling schemes
# using our method. This file is run on HPC cluster and the output .Rdata file is used for downstream analysis.

setwd("~/R/novogene_trimmed")
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.0",.libPaths()))
library("quantreg")
library("parallel")
source("src_funcs.r")
gene_keys <- read.csv(file = "Data_SampleKey_v3.csv", header = TRUE)
sample_index_18C <- which(gene_keys$Temperature == '18C')
time_18C <- (gene_keys$Day[sample_index_18C]-1)*24 + gene_keys$Hour[sample_index_18C]
sample_index_25C <- which(gene_keys$Temperature == '25C')
time_25C <- (gene_keys$Day[sample_index_25C]-1)*24 + gene_keys$Hour[sample_index_25C]

setwd("~/R/novogene_lambdalbd/00_SyntheticData")
SynData <- read.csv(file='Data_Synthetic_lambda0.50.csv',header=TRUE)
SynData[,2:5] <- log(SynData[,2:5]+1)

############################################################
############################################################
# Set parameters used in the fitting, 2 rep, 1hr
set.seed(0)
ntrials <- 1000
sig <- 0.05
fluc <- 0.01
nretry <- 5
fitted_time <- seq(-120,120,0.5)
time_comb <- as.vector(matrix(rep(SynData$time,2),nrow=2,byrow=TRUE))
n <- length(time_comb)
idx <- sapply(time_comb, function(x,y) which(y==x), y=SynData$time)
lambdasave_1hr <- matrix(NA,nrow=100,ncol=4)

for (k in seq(100)) {
  
  outdata <- matrix(data = NA, nrow = 4, ncol = 21)
  colnames(outdata) <- c('geneid','gene_name','int_25C','cos_25C','sin_25C','decaycos','decaysin',
                         'decayint','cos_18C','sin_18C','lambda','pval_cycle_25C','pval_cycle_18C',
                         'pval_cycle_diff','pval_phase_diff','pval_amp_diff','pval_median_diff',
                         'small_tpm','nboot','niter','convergence')
  outdata <- as.data.frame(outdata)
  outdata$geneid <- "NA"
  outdata$gene_name <- colnames(SynData)[2:5]
  outdata$small_tpm <- FALSE
  
  # Setup output matrix - AICc value in model selection
  outdata2 <- matrix(data = NA, nrow = 4, ncol = 25)
  colnames(outdata2) <- c('geneid','gene_name','type','AICc_01full','AICc_02full','AICc_07full',
                          'AICc_08full','AICc_11full','AICc_12full','AICc_13full','AICc_14full',
                          'AICc_15full','AICc_16full','AICc_01partial','AICc_02partial',
                          'AICc_07partial','AICc_08partial','AICc_11partial','AICc_12partial',
                          'AICc_13partial','AICc_14partial','AICc_15partial','AICc_16partial',
                          'model','TRUE_lambda')
  outdata2 <- as.data.frame(outdata2)
  outdata2$gene_name <- colnames(SynData)[2:5]
  
  # Setup output matrix - Parameter fit for the best model with minimum AICc
  outdata3 <- matrix(data = NA, nrow = 4, ncol = 11)
  colnames(outdata3) <- c('geneid','gene_name','V1','V2','V3','V4','V5','V6','V7','V8','V9')
  outdata3 <- as.data.frame(outdata3)
  outdata3$geneid <- "NA"
  outdata3$gene_name <- colnames(SynData)[2:5]
  
  for (i in seq(4)) {
    
    value <- SynData[idx,i+1]
    res <- sapply(value,function(x) rexp(1,rate=15/x)*sample(c(-1,1),size=1))
    
    temp_data <- data.frame(
      'value' = sapply(value+res,max,0),
      'time' = time_comb
    )
    
    #######
    # pre-process the data, generate an estimate of initial guess
    par_est <- get_guess(value = temp_data$value, time = temp_data$time,
                         coarse = seq(0.04,1,0.03), fine = seq(0.04,1,0.01))
    
    # nonlinear fit
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
    outdata$convergence[i] <- !errmsg
    outdata[i,3:11] <- nlrqmodel$par
    outdata$niter[i] <- as.vector(nlrqmodel$counts[1])
    
    #######
    # Wild Bootstrapping
    data_fitted <- data.frame(
      'value' = func(nlrqmodel$par,time_comb),
      'time' = time_comb
    )
    residual <- temp_data$value - data_fitted$value
    wbweights <- lapply(rep(length(time_comb),ntrials), get_weight)
    par_save_list <- mclapply(wbweights, wild_boot, fit_residual=residual, fit_data=data_fitted,
                              coarse=seq(0.04,1,0.03), fine=seq(0.04,1,0.01),
                              boot_nretry = nretry, boot_fluc = fluc,
                              mc.cores = 100)
    par_save <- do.call(rbind, par_save_list)
    outdata$nboot[i] <-  sum(!is.na(par_save[,1]))
    
    # significance
    cycle_para <- data.frame(
      x = par_save[,2],
      y = par_save[,3],
      p = par_save[,2] + par_save[,7],
      q = par_save[,3] + par_save[,8],
      a = par_save[,7],
      b = par_save[,8],
      f = sqrt((par_save[,2]+par_save[,7])^2+(par_save[,3]+par_save[,8])^2) -
        sqrt(par_save[,2]^2+par_save[,3]^2),
      g = par_save[,4] + par_save[,6] + par_save[,7]
    )
    cycle_para <- na.omit(cycle_para)
    ## cycle significance
    outdata$pval_cycle_25C[i] <- test_2D(Testdata=as.matrix(cycle_para[,1:2]))
    outdata$pval_cycle_18C[i] <- test_2D(Testdata=as.matrix(cycle_para[,3:4]))
    outdata$pval_cycle_diff[i] <- test_2D(Testdata=as.matrix(cycle_para[,5:6]))
    ## phase diff
    fitted_phase_diff <- get_angle(c(outdata[i,4],outdata[i,5],outdata[i,4]+outdata[i,9],outdata[i,5]+outdata[i,10]))
    outdata$pval_phase_diff[i] <- test_1D(Testdata=apply(as.matrix(cycle_para[,1:4]),1,get_angle))
    ## amp diff
    fitted_amp_diff <- sqrt((outdata[i,4]+outdata[i,9])^2+(outdata[i,5]+outdata[i,10])^2) -
      sqrt(outdata[i,4]^2+outdata[i,5]^2)
    outdata$pval_amp_diff[i] <- test_1D(Testdata=cycle_para[,7])
    ## median diff
    fitted_median_diff <- outdata[i,6]+outdata[i,8]+outdata[i,9]
    outdata$pval_median_diff[i] <- test_1D(Testdata=cycle_para[,8])
    
    ## full model v.s. partial model
    AICsave <- rep(NA,20)
    type_vec <- c('01f','02f','07f','08f','11f','12f','13f','14f','15f','16f',
                  '01p','02p','07p','08p','11p','12p','13p','14p','15p','16p')
    AICc_lambda <- rep(NA,length(type_vec))
    AICc_fitpar <- vector("list", 20)
    type_vec_num <- c(1,2,7,8,11,12,13,14,15)
    
    for (j in seq(length(type_vec_num))) {
      AICc_full <- get_fullmodel_AICc(data = temp_data, nretry = 5, type = type_vec_num[j])
      AICc_partial <- get_partialmodel_AICc(data = temp_data, type_vec_num[j])
      outdata2[i,j+3] <- AICc_full[[1]]
      AICc_lambda[j] <- AICc_full[[2]]
      outdata2[i,j+13] <-AICc_partial[[1]]
      AICc_lambda[j+10] <- AICc_partial[[2]]
      AICc_fitpar[[j]] <- AICc_full[[3]]
      AICc_fitpar[[j+10]] <- AICc_partial[[3]]
    }
    
    outdata2$AICc_16full[i] <- 2*9 - 2*(n*log(n/4/nlrqmodel$value)-n) + (2*9^2+2*9)/(n-9-1)
    AICc_lambda[10] <- nlrqmodel$par[9]
    AICc_fitpar[[10]] <- nlrqmodel$par
    
    AICc_16partial <- get_partialmodel_AICc(data = temp_data, type = 16)
    outdata2$AICc_16partial[i] <- AICc_16partial[[1]]
    AICc_lambda[20] <- AICc_16partial[[2]]
    AICc_fitpar[[20]] <- AICc_16partial[[3]]
    
    outdata2$type[i] <- (outdata$pval_cycle_25C[i]<=sig)*8 + (outdata$pval_cycle_18C[i]<=sig)*4 +
      (outdata$pval_cycle_diff[i]<=sig)*2 + (outdata$pval_median_diff[i]<=sig) + 1
    type_vec_ind <- which(outdata2[i,4:23] == min(outdata2[i,4:23]))
    outdata2$model[i] <- type_vec[type_vec_ind]
    outdata2$TRUE_lambda[i] <- AICc_lambda[type_vec_ind]
    
    outdata3[i,3:(2+length(AICc_fitpar[[type_vec_ind]]))] <- AICc_fitpar[[type_vec_ind]]
    
    message('1HR, ',i,'/4, ',k,'/100')
    
  }
  
  ind_concordant <- which(as.numeric(substr(outdata2$model,1,2))==outdata2$type &
                            apply(outdata2[,4:23],1,function(x) diff(sort(unlist(x)))[1]) >= 2)
  lambdasave_1hr[k,ind_concordant] <- outdata2$TRUE_lambda[ind_concordant]
  
}


############################################################
############################################################
# Set parameters used in the fitting, 2 rep, 2hr
set.seed(0)
ntrials <- 1000
sig <- 0.05
fluc <- 0.01
nretry <- 5
fitted_time <- seq(-120,120,0.5)
time_comb <- as.vector(matrix(rep(SynData$time[SynData$time%%2==0],2),nrow=2,byrow=TRUE))
n <- length(time_comb)
idx <- sapply(time_comb, function(x,y) which(y==x), y=SynData$time)
lambdasave_2hr <- matrix(NA,nrow=100,ncol=4)

for (k in seq(100)) {
  
  outdata <- matrix(data = NA, nrow = 4, ncol = 21)
  colnames(outdata) <- c('geneid','gene_name','int_25C','cos_25C','sin_25C','decaycos','decaysin',
                         'decayint','cos_18C','sin_18C','lambda','pval_cycle_25C','pval_cycle_18C',
                         'pval_cycle_diff','pval_phase_diff','pval_amp_diff','pval_median_diff',
                         'small_tpm','nboot','niter','convergence')
  outdata <- as.data.frame(outdata)
  outdata$geneid <- "NA"
  outdata$gene_name <- colnames(SynData)[2:5]
  outdata$small_tpm <- FALSE
  
  # Setup output matrix - AICc value in model selection
  outdata2 <- matrix(data = NA, nrow = 4, ncol = 25)
  colnames(outdata2) <- c('geneid','gene_name','type','AICc_01full','AICc_02full','AICc_07full',
                          'AICc_08full','AICc_11full','AICc_12full','AICc_13full','AICc_14full',
                          'AICc_15full','AICc_16full','AICc_01partial','AICc_02partial',
                          'AICc_07partial','AICc_08partial','AICc_11partial','AICc_12partial',
                          'AICc_13partial','AICc_14partial','AICc_15partial','AICc_16partial',
                          'model','TRUE_lambda')
  outdata2 <- as.data.frame(outdata2)
  outdata2$gene_name <- colnames(SynData)[2:5]
  
  # Setup output matrix - Parameter fit for the best model with minimum AICc
  outdata3 <- matrix(data = NA, nrow = 4, ncol = 11)
  colnames(outdata3) <- c('geneid','gene_name','V1','V2','V3','V4','V5','V6','V7','V8','V9')
  outdata3 <- as.data.frame(outdata3)
  outdata3$geneid <- "NA"
  outdata3$gene_name <- colnames(SynData)[2:5]
  
  for (i in seq(4)) {
    
    value <- SynData[idx,i+1]
    res <- sapply(value,function(x) rexp(1,rate=15/x)*sample(c(-1,1),size=1))
    
    temp_data <- data.frame(
      'value' = sapply(value+res,max,0),
      'time' = time_comb
    )
    
    #######
    # pre-process the data, generate an estimate of initial guess
    par_est <- get_guess(value = temp_data$value, time = temp_data$time,
                         coarse = seq(0.04,1,0.03), fine = seq(0.04,1,0.01))
    
    # nonlinear fit
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
    outdata$convergence[i] <- !errmsg
    outdata[i,3:11] <- nlrqmodel$par
    outdata$niter[i] <- as.vector(nlrqmodel$counts[1])
    
    #######
    # Wild Bootstrapping
    data_fitted <- data.frame(
      'value' = func(nlrqmodel$par,time_comb),
      'time' = time_comb
    )
    residual <- temp_data$value - data_fitted$value
    wbweights <- lapply(rep(length(time_comb),ntrials), get_weight)
    par_save_list <- mclapply(wbweights, wild_boot, fit_residual=residual, fit_data=data_fitted,
                              coarse=seq(0.04,1,0.03), fine=seq(0.04,1,0.01),
                              boot_nretry = nretry, boot_fluc = fluc,
                              mc.cores = 100)
    par_save <- do.call(rbind, par_save_list)
    outdata$nboot[i] <-  sum(!is.na(par_save[,1]))
    
    # significance
    cycle_para <- data.frame(
      x = par_save[,2],
      y = par_save[,3],
      p = par_save[,2] + par_save[,7],
      q = par_save[,3] + par_save[,8],
      a = par_save[,7],
      b = par_save[,8],
      f = sqrt((par_save[,2]+par_save[,7])^2+(par_save[,3]+par_save[,8])^2) -
        sqrt(par_save[,2]^2+par_save[,3]^2),
      g = par_save[,4] + par_save[,6] + par_save[,7]
    )
    cycle_para <- na.omit(cycle_para)
    ## cycle significance
    outdata$pval_cycle_25C[i] <- test_2D(Testdata=as.matrix(cycle_para[,1:2]))
    outdata$pval_cycle_18C[i] <- test_2D(Testdata=as.matrix(cycle_para[,3:4]))
    outdata$pval_cycle_diff[i] <- test_2D(Testdata=as.matrix(cycle_para[,5:6]))
    ## phase diff
    fitted_phase_diff <- get_angle(c(outdata[i,4],outdata[i,5],outdata[i,4]+outdata[i,9],outdata[i,5]+outdata[i,10]))
    outdata$pval_phase_diff[i] <- test_1D(Testdata=apply(as.matrix(cycle_para[,1:4]),1,get_angle))
    ## amp diff
    fitted_amp_diff <- sqrt((outdata[i,4]+outdata[i,9])^2+(outdata[i,5]+outdata[i,10])^2) -
      sqrt(outdata[i,4]^2+outdata[i,5]^2)
    outdata$pval_amp_diff[i] <- test_1D(Testdata=cycle_para[,7])
    ## median diff
    fitted_median_diff <- outdata[i,6]+outdata[i,8]+outdata[i,9]
    outdata$pval_median_diff[i] <- test_1D(Testdata=cycle_para[,8])
    
    ## full model v.s. partial model
    AICsave <- rep(NA,20)
    type_vec <- c('01f','02f','07f','08f','11f','12f','13f','14f','15f','16f',
                  '01p','02p','07p','08p','11p','12p','13p','14p','15p','16p')
    AICc_lambda <- rep(NA,length(type_vec))
    AICc_fitpar <- vector("list", 20)
    type_vec_num <- c(1,2,7,8,11,12,13,14,15)
    
    for (j in seq(length(type_vec_num))) {
      AICc_full <- get_fullmodel_AICc(data = temp_data, nretry = 5, type = type_vec_num[j])
      AICc_partial <- get_partialmodel_AICc(data = temp_data, type_vec_num[j])
      outdata2[i,j+3] <- AICc_full[[1]]
      AICc_lambda[j] <- AICc_full[[2]]
      outdata2[i,j+13] <-AICc_partial[[1]]
      AICc_lambda[j+10] <- AICc_partial[[2]]
      AICc_fitpar[[j]] <- AICc_full[[3]]
      AICc_fitpar[[j+10]] <- AICc_partial[[3]]
    }
    
    outdata2$AICc_16full[i] <- 2*9 - 2*(n*log(n/4/nlrqmodel$value)-n) + (2*9^2+2*9)/(n-9-1)
    AICc_lambda[10] <- nlrqmodel$par[9]
    AICc_fitpar[[10]] <- nlrqmodel$par
    
    AICc_16partial <- get_partialmodel_AICc(data = temp_data, type = 16)
    outdata2$AICc_16partial[i] <- AICc_16partial[[1]]
    AICc_lambda[20] <- AICc_16partial[[2]]
    AICc_fitpar[[20]] <- AICc_16partial[[3]]
    
    outdata2$type[i] <- (outdata$pval_cycle_25C[i]<=sig)*8 + (outdata$pval_cycle_18C[i]<=sig)*4 +
      (outdata$pval_cycle_diff[i]<=sig)*2 + (outdata$pval_median_diff[i]<=sig) + 1
    type_vec_ind <- which(outdata2[i,4:23] == min(outdata2[i,4:23]))
    outdata2$model[i] <- type_vec[type_vec_ind]
    outdata2$TRUE_lambda[i] <- AICc_lambda[type_vec_ind]
    
    outdata3[i,3:(2+length(AICc_fitpar[[type_vec_ind]]))] <- AICc_fitpar[[type_vec_ind]]
    
    message('2HR, ',i,'/4, ',k,'/100')
    
  }
  
  ind_concordant <- which(as.numeric(substr(outdata2$model,1,2))==outdata2$type &
                            apply(outdata2[,4:23],1,function(x) diff(sort(unlist(x)))[1]) >= 2)
  lambdasave_2hr[k,ind_concordant] <- outdata2$TRUE_lambda[ind_concordant]
  
}


############################################################
############################################################
# Set parameters used in the fitting, our sampling scheme
set.seed(0)
ntrials <- 1000
sig <- 0.05
fluc <- 0.01
nretry <- 5
fitted_time <- seq(-120,120,0.5)
time_comb <- c(time_25C-120,time_18C)
idx <- sapply(time_comb, function(x,y) which(y==x), y=SynData$time)
n <- length(time_comb)
lambdasave_exp <- matrix(NA,nrow=100,ncol=4)

for (k in seq(100)) {
  
  outdata <- matrix(data = NA, nrow = 4, ncol = 21)
  colnames(outdata) <- c('geneid','gene_name','int_25C','cos_25C','sin_25C','decaycos','decaysin',
                         'decayint','cos_18C','sin_18C','lambda','pval_cycle_25C','pval_cycle_18C',
                         'pval_cycle_diff','pval_phase_diff','pval_amp_diff','pval_median_diff',
                         'small_tpm','nboot','niter','convergence')
  outdata <- as.data.frame(outdata)
  outdata$geneid <- "NA"
  outdata$gene_name <- colnames(SynData)[2:5]
  outdata$small_tpm <- FALSE
  
  # Setup output matrix - AICc value in model selection
  outdata2 <- matrix(data = NA, nrow = 4, ncol = 25)
  colnames(outdata2) <- c('geneid','gene_name','type','AICc_01full','AICc_02full','AICc_07full',
                          'AICc_08full','AICc_11full','AICc_12full','AICc_13full','AICc_14full',
                          'AICc_15full','AICc_16full','AICc_01partial','AICc_02partial',
                          'AICc_07partial','AICc_08partial','AICc_11partial','AICc_12partial',
                          'AICc_13partial','AICc_14partial','AICc_15partial','AICc_16partial',
                          'model','TRUE_lambda')
  outdata2 <- as.data.frame(outdata2)
  outdata2$gene_name <- colnames(SynData)[2:5]
  
  # Setup output matrix - Parameter fit for the best model with minimum AICc
  outdata3 <- matrix(data = NA, nrow = 4, ncol = 11)
  colnames(outdata3) <- c('geneid','gene_name','V1','V2','V3','V4','V5','V6','V7','V8','V9')
  outdata3 <- as.data.frame(outdata3)
  outdata3$geneid <- "NA"
  outdata3$gene_name <- colnames(SynData)[2:5]
  
  for (i in seq(4)) {
    
    value <- SynData[idx,i+1]
    res <- sapply(value,function(x) rexp(1,rate=15/x)*sample(c(-1,1),size=1))
    
    temp_data <- data.frame(
      'value' = sapply(value+res,max,0),
      'time' = time_comb
    )
    
    #######
    # pre-process the data, generate an estimate of initial guess
    par_est <- get_guess(value = temp_data$value, time = temp_data$time,
                         coarse = seq(0.04,1,0.03), fine = seq(0.04,1,0.01))
    
    # nonlinear fit
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
    outdata$convergence[i] <- !errmsg
    outdata[i,3:11] <- nlrqmodel$par
    outdata$niter[i] <- as.vector(nlrqmodel$counts[1])
    
    #######
    # Wild Bootstrapping
    data_fitted <- data.frame(
      'value' = func(nlrqmodel$par,time_comb),
      'time' = time_comb
    )
    residual <- temp_data$value - data_fitted$value
    wbweights <- lapply(rep(length(time_comb),ntrials), get_weight)
    par_save_list <- mclapply(wbweights, wild_boot, fit_residual=residual, fit_data=data_fitted,
                              coarse=seq(0.04,1,0.03), fine=seq(0.04,1,0.01),
                              boot_nretry = nretry, boot_fluc = fluc,
                              mc.cores = 100)
    par_save <- do.call(rbind, par_save_list)
    outdata$nboot[i] <-  sum(!is.na(par_save[,1]))
    
    # significance
    cycle_para <- data.frame(
      x = par_save[,2],
      y = par_save[,3],
      p = par_save[,2] + par_save[,7],
      q = par_save[,3] + par_save[,8],
      a = par_save[,7],
      b = par_save[,8],
      f = sqrt((par_save[,2]+par_save[,7])^2+(par_save[,3]+par_save[,8])^2) -
        sqrt(par_save[,2]^2+par_save[,3]^2),
      g = par_save[,4] + par_save[,6] + par_save[,7]
    )
    cycle_para <- na.omit(cycle_para)
    ## cycle significance
    outdata$pval_cycle_25C[i] <- test_2D(Testdata=as.matrix(cycle_para[,1:2]))
    outdata$pval_cycle_18C[i] <- test_2D(Testdata=as.matrix(cycle_para[,3:4]))
    outdata$pval_cycle_diff[i] <- test_2D(Testdata=as.matrix(cycle_para[,5:6]))
    ## phase diff
    fitted_phase_diff <- get_angle(c(outdata[i,4],outdata[i,5],outdata[i,4]+outdata[i,9],outdata[i,5]+outdata[i,10]))
    outdata$pval_phase_diff[i] <- test_1D(Testdata=apply(as.matrix(cycle_para[,1:4]),1,get_angle))
    ## amp diff
    fitted_amp_diff <- sqrt((outdata[i,4]+outdata[i,9])^2+(outdata[i,5]+outdata[i,10])^2) -
      sqrt(outdata[i,4]^2+outdata[i,5]^2)
    outdata$pval_amp_diff[i] <- test_1D(Testdata=cycle_para[,7])
    ## median diff
    fitted_median_diff <- outdata[i,6]+outdata[i,8]+outdata[i,9]
    outdata$pval_median_diff[i] <- test_1D(Testdata=cycle_para[,8])
    
    ## full model v.s. partial model
    AICsave <- rep(NA,20)
    type_vec <- c('01f','02f','07f','08f','11f','12f','13f','14f','15f','16f',
                  '01p','02p','07p','08p','11p','12p','13p','14p','15p','16p')
    AICc_lambda <- rep(NA,length(type_vec))
    AICc_fitpar <- vector("list", 20)
    type_vec_num <- c(1,2,7,8,11,12,13,14,15)
    
    for (j in seq(length(type_vec_num))) {
      AICc_full <- get_fullmodel_AICc(data = temp_data, nretry = 5, type = type_vec_num[j])
      AICc_partial <- get_partialmodel_AICc(data = temp_data, type_vec_num[j])
      outdata2[i,j+3] <- AICc_full[[1]]
      AICc_lambda[j] <- AICc_full[[2]]
      outdata2[i,j+13] <-AICc_partial[[1]]
      AICc_lambda[j+10] <- AICc_partial[[2]]
      AICc_fitpar[[j]] <- AICc_full[[3]]
      AICc_fitpar[[j+10]] <- AICc_partial[[3]]
    }
    
    outdata2$AICc_16full[i] <- 2*9 - 2*(n*log(n/4/nlrqmodel$value)-n) + (2*9^2+2*9)/(n-9-1)
    AICc_lambda[10] <- nlrqmodel$par[9]
    AICc_fitpar[[10]] <- nlrqmodel$par
    
    AICc_16partial <- get_partialmodel_AICc(data = temp_data, type = 16)
    outdata2$AICc_16partial[i] <- AICc_16partial[[1]]
    AICc_lambda[20] <- AICc_16partial[[2]]
    AICc_fitpar[[20]] <- AICc_16partial[[3]]
    
    outdata2$type[i] <- (outdata$pval_cycle_25C[i]<=sig)*8 + (outdata$pval_cycle_18C[i]<=sig)*4 +
      (outdata$pval_cycle_diff[i]<=sig)*2 + (outdata$pval_median_diff[i]<=sig) + 1
    type_vec_ind <- which(outdata2[i,4:23] == min(outdata2[i,4:23]))
    outdata2$model[i] <- type_vec[type_vec_ind]
    outdata2$TRUE_lambda[i] <- AICc_lambda[type_vec_ind]
    
    outdata3[i,3:(2+length(AICc_fitpar[[type_vec_ind]]))] <- AICc_fitpar[[type_vec_ind]]
    
    message('EXP, ',i,'/4, ',k,'/100')
    
  }
  
  ind_concordant <- which(as.numeric(substr(outdata2$model,1,2))==outdata2$type &
                            apply(outdata2[,4:23],1,function(x) diff(sort(unlist(x)))[1]) >= 2)
  lambdasave_exp[k,ind_concordant] <- outdata2$TRUE_lambda[ind_concordant]
  
}

save.image(file='00.03SyntheticData_lambda0.50.Rdata')

