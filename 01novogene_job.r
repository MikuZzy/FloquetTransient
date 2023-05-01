# This is the main script that performs our method on the experimental data. This file is run on HPC
# cluster (mclapply function) and the output .Rdata and .csv files are used for downstream analysis.

setwd("~/R/novogene_trimmed")
library("quantreg")
library("parallel")
source("src_funcs.r")

###########################################################
############################################################
# Get Data
Flybase_IDs_CirRhy <- read.table("Data_FlyBase_IDs_CirRhy.txt")
Flybase_IDs_Resp2Temp <- read.table("Data_FlyBase_IDs_Resp2Temp.txt")
Flybase_IDs_TransReg <- read.table("Data_FlyBase_IDs_TransReg.txt")
novogene2021 <- read.table(file = "Data_novogene2021_trimmed_dedup_tpm.tsv", header = TRUE)
gene_keys <- read.csv(file = "Data_SampleKey_v3.csv", header = TRUE)

# 18C Data from the current experiment
sample_index_18C <- which(gene_keys$Temperature == '18C')
time_18C <- (gene_keys$Day[sample_index_18C]-1)*24 + gene_keys$Hour[sample_index_18C]
allgenes_data_18C <- novogene2021[,sample_index_18C+2]
ngene <- nrow(allgenes_data_18C)
mindetect_18C <- rep(NA, ncol(allgenes_data_18C))
for (i in seq(ncol(allgenes_data_18C))) {
  mindetect_18C[i] <- min(allgenes_data_18C[which(allgenes_data_18C[,i]!=0),i])
}

# 25C Data from the current experiment
sample_index_25C <- which(gene_keys$Temperature == '25C')
time_25C <- (gene_keys$Day[sample_index_25C]-1)*24 + gene_keys$Hour[sample_index_25C]
allgenes_data_25C <- novogene2021[,sample_index_25C+2]
mindetect_25C <- rep(NA, ncol(allgenes_data_25C))
for (i in seq(ncol(allgenes_data_25C))) {
  mindetect_25C[i] <- min(allgenes_data_25C[which(allgenes_data_25C[,i]!=0),i])
}

# ZT for all data
time_comb <- c(time_25C-120,time_18C)
n <- length(time_comb)

############################################################
############################################################
# Setup output matrix - parameter fit for the full model
outdata <- matrix(data = NA, nrow = ngene, ncol = 21)
colnames(outdata) <- c('geneid','gene_name','int_25C','cos_25C','sin_25C','decaycos','decaysin',
                       'decayint','cos_18C','sin_18C','lambda','pval_cycle_25C','pval_cycle_18C',
                       'pval_cycle_diff','pval_phase_diff','pval_amp_diff','pval_median_diff',
                       'small_tpm','nboot','niter','convergence')
outdata <- as.data.frame(outdata)
outdata$geneid <- novogene2021$geneid
outdata$gene_name <- novogene2021$gene_name
outdata$small_tpm <- FALSE

# Setup output matrix - AICc value in model selection
outdata2 <- matrix(data = NA, nrow = ngene, ncol = 25)
colnames(outdata2) <- c('geneid','gene_name','type','AICc_01full','AICc_02full','AICc_07full',
                        'AICc_08full','AICc_11full','AICc_12full','AICc_13full','AICc_14full',
                        'AICc_15full','AICc_16full','AICc_01partial','AICc_02partial',
                        'AICc_07partial','AICc_08partial','AICc_11partial','AICc_12partial',
                        'AICc_13partial','AICc_14partial','AICc_15partial','AICc_16partial',
                        'model','TRUE_lambda')
outdata2 <- as.data.frame(outdata2)
outdata2$geneid <- novogene2021$geneid
outdata2$gene_name <- novogene2021$gene_name

# Setup output matrix - Parameter fit for the best model with minimum AICc
outdata3 <- matrix(data = NA, nrow = ngene, ncol = 11)
colnames(outdata3) <- c('geneid','gene_name','V1','V2','V3','V4','V5','V6','V7','V8','V9')
outdata3 <- as.data.frame(outdata3)
outdata3$geneid <- novogene2021$geneid
outdata3$gene_name <- novogene2021$gene_name

############################################################
############################################################
# Set parameters used in the fitting
set.seed(0)
ntrials <- 1000
sig <- 0.05
fluc <- 0.01
nretry <- 5
fitted_time <- seq(-120,120,0.5)

# Main loop - fit
for (i in seq(ngene)) {
  
  value_18C <- unlist(allgenes_data_18C[i,])
  value_25C <- unlist(allgenes_data_25C[i,])
  if (median(value_18C) <= 5 & median(value_25C) <= 5) {
    outdata$small_tpm[i] <- TRUE
  } else {
    value_18C[which(value_18C==0)] <- mindetect_18C[which(value_18C==0)]
    value_25C[which(value_25C==0)] <- mindetect_25C[which(value_25C==0)]
    temp_data <- data.frame(
      'value' = c(value_25C,value_18C),
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

  }
  
  message(i)
  
}

############################################################
############################################################
# Save outputs
setwd("~/R/novogene_trimmed")
write.csv(outdata, file = '01Result_fullparafit.csv')
write.csv(outdata2, file = '01Result_modelAICc.csv')
write.csv(outdata3, file = '01Result_bestmodel.csv')

save.image(file = '01novogene_job.Rdata')

############################################################
############################################################
# Filter genes with small TPM
ind_notNA <- which(outdata$small_tpm == FALSE)
outdata_filtered <- outdata[ind_notNA,]
outdata2_filtered <- outdata2[ind_notNA,]
outdata3_filtered <- outdata3[ind_notNA,]
write.csv(outdata_filtered,file='01Result_filtered_fullparafit.csv')
write.csv(outdata2_filtered,file='01Result_filtered_modelAICc.csv')
write.csv(outdata3_filtered,file='01Result_filtered_bestmodel.csv')

############################################################
############################################################
# Sort out genes with concordant estimated types
colvec <- c('01f','01p','02f','02p','07f','07p','08f','08p','11f','11p','12f','12p',
            '13f','13p','14f','14p','15f','15p','16f','16p')
rowvec <- seq(16)
clf <- matrix(0,nrow=16,ncol=20)
colnames(clf) <- colvec
rownames(clf) <- rowvec
ksave <- c()
count <- 0
for (k in seq(nrow(outdata_filtered))) {
  i <- which(rowvec == outdata2_filtered$type[k])
  j <- which(colvec == outdata2_filtered$model[k])
  AICc_vec <- as.vector(unlist(outdata2_filtered[k,4:23]))
  AICc_vec <- AICc_vec[order(AICc_vec)]
  dAICc <- diff(AICc_vec)[1]
  clf[i,j] <- clf[i,j] + 1
  if (as.numeric(substr(outdata2_filtered$model[k],1,2))==outdata2_filtered$type[k]) {
    count <- count + 1
    if (dAICc>2) {
       ksave <- append(ksave,k)
    }
  }
  message(k,'/',nrow(outdata_filtered))
}
write.csv(clf,file='01Result_concordant_overview.csv')
outdata2_sel_cont <- outdata2_filtered[ksave,]
write.csv(outdata2_sel_cont,file='01Result_concordant_details.csv')

############################################################
############################################################
# Sort out reproducible genes
setwd("~/R/novogene_trimmed/00_Reproducibility")
result_filter <- read.csv(file="00ReproducibleGenes.csv")
setwd("~/R/novogene_trimmed")
outdata2_sel_cont_filter <- outdata2_sel_cont[which(outdata2_sel_cont$geneid %in% result_filter$geneid),]
write.csv(outdata2_sel_cont_filter,file='01Result_concordant_details_repro.csv')

outdata_repro <- outdata_filtered[which(outdata_filtered$geneid %in% result_filter$geneid),]
outdata2_repro <- outdata2_filtered[which(outdata2_filtered$geneid %in% result_filter$geneid),]

############################################################
############################################################
# Sort out genes from three GO categories
CirRhy <- results_from_FlybaseIDs(IDList=Flybase_IDs_CirRhy$V1,
                                  result1=outdata_repro,result2=outdata2_repro)
write.csv(CirRhy,file="01Result_GOCirRhy.csv")

Resp2Temp <- results_from_FlybaseIDs(IDList=Flybase_IDs_Resp2Temp$V1,
                                     result1=outdata_repro,result2=outdata2_repro)
write.csv(Resp2Temp,file="01Result_GOResp2Temp.csv")

TransReg <- results_from_FlybaseIDs(IDList=Flybase_IDs_TransReg$V1,
                                    result1=outdata_repro,result2=outdata2_repro)
write.csv(TransReg,file="01Result_GOTransReg.csv")


save.image(file = '01novogene_job.Rdata')
