# This file creates synthetic data based on the full model fit for gene per and apply our method to them. 
# Residual distributions' MAD are varied and for each MAD there are n = 1000 trials. This file is run on HPC
# cluster (mclapply function) and the output .Rdata and .csv files are used for downstream analysis.

setwd("~/R/novogene_trimmed")
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.0",.libPaths()))
library('quantreg')
library('parallel')
source("src_funcs.r")

tpm_data <- read.table(file = 'Data_novogene2021_trimmed_dedup_tpm.tsv', header = TRUE)
gene_keys <- read.csv(file = 'Data_SampleKey_v3.csv', header = TRUE)
results_01 <- read.csv(file='01Result_fullparafit.csv', header = TRUE)
results_02 <- read.csv(file='01Result_modelAICc.csv', header = TRUE)
results_03 <- read.csv(file='01Result_bestmodel.csv', header = TRUE)

sample_index_18C <- which(gene_keys$Temperature == '18C')
time_18C <- (gene_keys$Day[sample_index_18C]-1)*24 + gene_keys$Hour[sample_index_18C]
sample_index_25C <- which(gene_keys$Temperature == '25C')
time_25C <- (gene_keys$Day[sample_index_25C]-1)*24 + gene_keys$Hour[sample_index_25C]
time_comb <- c(time_25C-120,time_18C)
n <- length(time_comb)

set.seed(0)
ntrials <- 1000
sig <- 0.05
fluc <- 0.01
nretry <- 5
fitted_time <- seq(-120,120,0.5)

p <- 490
par_value <- as.vector(as.numeric(results_01[p,4:12]))
fitted <- func_fullmodel_type(par_value,time_comb,type=16)
rawtpm <- as.vector(c(unlist(tpm_data[p,sample_index_25C+2]),unlist(tpm_data[p,sample_index_18C+2])))
MAD <- mean(abs(rawtpm-fitted))
b_value <- seq(0.5,5,0.5)

FNs_nosel <- matrix(data=NA, nrow=10, ncol=16)
colnames(FNs_nosel) <- c('Type01','Type02','Type03','Type04','Type05','Type06','Type07','Type08',
                         'Type09','Type10','Type11','Type12','Type13','Type14','Type15','Type16')
FNs_nosel <- as.data.frame(FNs_nosel)

FNs_withsel <- matrix(data=NA, nrow=10, ncol=21)
colnames(FNs_withsel) <- c('Type01f','Type02f','Type01p','Type02p','Type07f','Type08f','Type07p','Type08p',
                           'Type11f','Type12f','Type11p','Type12p','Type13f','Type14f','Type15f','Type16f',
                           'Type13p','Type14p','Type15p','Type16p','other')
FNs_withsel <- as.data.frame(FNs_withsel)

for (k in seq(10)) {
  
  b <- b_value[k]
  outdata <- matrix(data = NA, nrow = 1000, ncol = 21)
  colnames(outdata) <- c('geneid','gene_name','int_25C','cos_25C','sin_25C','decaycos','decaysin',
                         'decayint','cos_18C','sin_18C','lambda','pval_cycle_25C','pval_cycle_18C',
                         'pval_cycle_diff','pval_phase_diff','pval_amp_diff','pval_median_diff',
                         'small_tpm','nboot','niter','convergence')
  outdata <- as.data.frame(outdata)
  
  outdata2 <- matrix(data = NA, nrow = 1000, ncol = 26)
  colnames(outdata2) <- c('geneid','gene_name','type','AICc_01full','AICc_02full','AICc_07full',
                          'AICc_08full','AICc_11full','AICc_12full','AICc_13full','AICc_14full',
                          'AICc_15full','AICc_16full','AICc_01partial','AICc_02partial',
                          'AICc_07partial','AICc_08partial','AICc_11partial','AICc_12partial',
                          'AICc_13partial','AICc_14partial','AICc_15partial','AICc_16partial',
                          'model','TRUE_lambda','dAICc')
  outdata2 <- as.data.frame(outdata2)
  
  for (i in seq(1000)) {
    res_value <- rexp(n,rate=1/b) * sample(c(-1,1),size=n,replace=TRUE)
    temp_data <- data.frame(
      'value' = fitted + res_value,
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
                              mc.cores = 100L)
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
    outdata2$dAICc[i] <- diff(sort(unlist(outdata2[i,4:23])))[1]
    
    message(k,' ',i,'/1000')
  }
  
  FNs_nosel$Type01[k] <- sum(outdata2$type == 1, na.rm=TRUE)
  FNs_nosel$Type02[k] <- sum(outdata2$type == 2, na.rm=TRUE)
  FNs_nosel$Type03[k] <- sum(outdata2$type == 3, na.rm=TRUE)
  FNs_nosel$Type04[k] <- sum(outdata2$type == 4, na.rm=TRUE)
  FNs_nosel$Type05[k] <- sum(outdata2$type == 5, na.rm=TRUE)
  FNs_nosel$Type06[k] <- sum(outdata2$type == 6, na.rm=TRUE)
  FNs_nosel$Type07[k] <- sum(outdata2$type == 7, na.rm=TRUE)
  FNs_nosel$Type08[k] <- sum(outdata2$type == 8, na.rm=TRUE)
  FNs_nosel$Type09[k] <- sum(outdata2$type == 9, na.rm=TRUE)
  FNs_nosel$Type10[k] <- sum(outdata2$type == 10, na.rm=TRUE)
  FNs_nosel$Type11[k] <- sum(outdata2$type == 11, na.rm=TRUE)
  FNs_nosel$Type12[k] <- sum(outdata2$type == 12, na.rm=TRUE)
  FNs_nosel$Type13[k] <- sum(outdata2$type == 13, na.rm=TRUE)
  FNs_nosel$Type14[k] <- sum(outdata2$type == 14, na.rm=TRUE)
  FNs_nosel$Type15[k] <- sum(outdata2$type == 15, na.rm=TRUE)
  FNs_nosel$Type16[k] <- sum(outdata2$type == 16, na.rm=TRUE)
  
  FNs_withsel$Type01f[k] <- sum(outdata2$type == 1 & outdata2$model == '01f' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type02f[k] <- sum(outdata2$type == 2 & outdata2$model == '02f' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type07f[k] <- sum(outdata2$type == 7 & outdata2$model == '07f' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type08f[k] <- sum(outdata2$type == 8 & outdata2$model == '08f' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type11f[k] <- sum(outdata2$type == 11 & outdata2$model == '11f' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type12f[k] <- sum(outdata2$type == 12 & outdata2$model == '12f' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type13f[k] <- sum(outdata2$type == 13 & outdata2$model == '13f' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type14f[k] <- sum(outdata2$type == 14 & outdata2$model == '14f' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type15f[k] <- sum(outdata2$type == 15 & outdata2$model == '15f' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type16f[k] <- sum(outdata2$type == 16 & outdata2$model == '16f' & outdata2$dAICc >= 2, na.rm=TRUE)
  
  FNs_withsel$Type01p[k] <- sum(outdata2$type == 1 & outdata2$model == '01p' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type02p[k] <- sum(outdata2$type == 2 & outdata2$model == '02p' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type07p[k] <- sum(outdata2$type == 7 & outdata2$model == '07p' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type08p[k] <- sum(outdata2$type == 8 & outdata2$model == '08p' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type11p[k] <- sum(outdata2$type == 11 & outdata2$model == '11p' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type12p[k] <- sum(outdata2$type == 12 & outdata2$model == '12p' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type13p[k] <- sum(outdata2$type == 13 & outdata2$model == '13p' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type14p[k] <- sum(outdata2$type == 14 & outdata2$model == '14p' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type15p[k] <- sum(outdata2$type == 15 & outdata2$model == '15p' & outdata2$dAICc >= 2, na.rm=TRUE)
  FNs_withsel$Type16p[k] <- sum(outdata2$type == 16 & outdata2$model == '16p' & outdata2$dAICc >= 2, na.rm=TRUE)
  
  FNs_withsel$other[k] <- 1000 - sum(FNs_withsel[k,1:20])
}

setwd("~/R/novogene_trimmed/02_SimulatedFNs")
write.csv(FNs_nosel, file = '02Result_PER_nosel.csv')
write.csv(FNs_withsel, file = '02Result_PER_withsel.csv')
save.image(file = '02FalseNegatives_PER.Rdata')
