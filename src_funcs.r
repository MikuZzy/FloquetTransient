# This is the main source function script, and most scripts need to import this source file.

####################
# full piece-wise nonlinear model, used in nlrq
func_new <- function(intercept,cos_25C,sin_25C,decaycos,decaysin,decayint,cos_18C,sin_18C,lambda,time) {
  ind1 <- which(time < 0)
  ind2 <- which(time >= 0)
  result <- time * 0
  result[ind1] <- intercept +
    cos_25C * cos(time[ind1]*pi/12) +
    sin_25C * sin(time[ind1]*pi/12)
  result[ind2] <- intercept +
    cos_25C * cos(time[ind2]*pi/12) +
    sin_25C * sin(time[ind2]*pi/12) +
    decaycos * (exp(-lambda*time[ind2]) * cos(time[ind2]*pi/12) - 1) +
    decaysin * (exp(-lambda*time[ind2]) * sin(time[ind2]*pi/12) - 0) +
    decayint * (exp(-lambda*time[ind2]) - 1) +
    cos_18C * (cos(time[ind2]*pi/12) - 1) +
    sin_18C * (sin(time[ind2]*pi/12) - 0)
  return(result)
}

# full piece-wise nonlinear model, used in plotting and optim
func <- function(par, time) {
  ind1 <- which(time < 0)
  ind2 <- which(time >= 0)
  result <- time * 0
  result[ind1] <- par[1] +
    par[2] * cos(time[ind1]*pi/12) +
    par[3] * sin(time[ind1]*pi/12)
  result[ind2] <- par[1] +
    par[2] * cos(time[ind2]*pi/12) +
    par[3] * sin(time[ind2]*pi/12) +
    par[4] * (exp(-par[9]*time[ind2]) * cos(time[ind2]*pi/12) - 1) +
    par[5] * (exp(-par[9]*time[ind2]) * sin(time[ind2]*pi/12) - 0) +
    par[6] * (exp(-par[9]*time[ind2]) - 1) +
    par[7] * (cos(time[ind2]*pi/12) - 1) +
    par[8] * (sin(time[ind2]*pi/12) - 0)
  return(result)
}

# objective function to minimize the full model, used in optim
obj <- function(par, value, time) {
  AD <- value - func(par, time)
  result <- 1/2 * sum(AD*sign(AD))
  return(result)
}

# gradient of full model's objective function, used in optim
grad <- function(par, value, time) {
  abs_fac <- sign(value - func(par, time))
  result <- rep(NA, 9)
  result[1] <- -1/2 * sum(abs_fac * 1)
  result[2] <- -1/2 * sum(abs_fac * cos(time*pi/12))
  result[3] <- -1/2 * sum(abs_fac * sin(time*pi/12))
  result[4] <- -1/2 * sum(abs_fac * (exp(-par[9]*time)*cos(time*pi/12)-1)*(time>=0))
  result[5] <- -1/2 * sum(abs_fac * (exp(-par[9]*time)*sin(time*pi/12)-0)*(time>=0))
  result[6] <- -1/2 * sum(abs_fac * (exp(-par[9]*time)-1)*(time>=0))
  result[7] <- -1/2 * sum(abs_fac * (cos(time*pi/12)-1)*(time>=0))
  result[8] <- -1/2 * sum(abs_fac * (sin(time*pi/12)-0)*(time>=0))
  result[9] <- -1/2 * sum(abs_fac * -time*exp(-par[9]*time)*(par[4]*cos(time*pi/12)+par[5]*sin(time*pi/12)+par[6])
                          *(time>=0))
  return(result)
}

# partial piece-wise nonlinear model, used in plotting
func_partial <- function(par, time) {
  ind1 <- which(time < 0)
  ind2 <- which(time >= 0)
  result <- time * 0
  result[ind1] <- par[1] +
    par[2] * cos(time[ind1]*pi/12) +
    par[3] * sin(time[ind1]*pi/12)
  result[ind2] <- par[1] +
    par[2] * cos(time[ind2]*pi/12) +
    par[3] * sin(time[ind2]*pi/12) +
    par[4] * (cos(time[ind2]*pi/12) - 1) +
    par[5] * (sin(time[ind2]*pi/12) - 0)
  return(result)
}

####################
# function to generate an initial parameter guess using sweep
get_guess <- function(value, time, coarse = seq(0.04,1,0.03), fine = NULL) {
  lambda_range <- coarse
  objective_storage <- 0*lambda_range
  objective_save <- Inf
  ind_save <- 0
  model_save <- NA
  for (j in seq(length(lambda_range))) {
    lambda <- lambda_range[j]
    temp_data <- data.frame(
      'value' = value,
      'time' = time,
      'cos_25C' = cos(time*pi/12),
      'sin_25C' = sin(time*pi/12),
      'decaycos' = (exp(-lambda*time) * cos(time*pi/12) - 1) * (time >= 0),
      'decaysin' = (exp(-lambda*time) * sin(time*pi/12) - 0) * (time >= 0),
      'decayint' = (exp(-lambda*time) - 1) * (time >= 0),
      'cos_18C' = (cos(time*pi/12) - 1) * (time >= 0),
      'sin_18C' = (sin(time*pi/12) - 0) * (time >= 0)
    )
    rqmodel <- rq(value ~ cos_25C + sin_25C + decaycos + decaysin + decayint + cos_18C + sin_18C + 1,
                  data = temp_data, tau = 0.5)
    objective <- rqmodel$rho
    objective_storage[j] <- objective
    if (objective < objective_save) {
      ind_save <- j
      objective_save <- objective
      model_save <- rqmodel
    }
  }
  if (!is.null(fine) & (lambda_range[ind_save] == max(lambda_range) || lambda_range[ind_save] == min(lambda_range))) {
    lambda_range <- fine
    objective_storage <- 0*lambda_range
    objective_save <- Inf
    ind_save <- 0
    model_save <- NA
    for (j in seq(length(lambda_range))) {
      lambda <- lambda_range[j]
      temp_data <- data.frame(
        'value' = value,
        'time' = time,
        'cos_25C' = cos(time*pi/12),
        'sin_25C' = sin(time*pi/12),
        'decaycos' = (exp(-lambda*time) * cos(time*pi/12) - 1) * (time >= 0),
        'decaysin' = (exp(-lambda*time) * sin(time*pi/12) - 0) * (time >= 0),
        'decayint' = (exp(-lambda*time) - 1) * (time >= 0),
        'cos_18C' = (cos(time*pi/12) - 1) * (time >= 0),
        'sin_18C' = (sin(time*pi/12) - 0) * (time >= 0)
      )
      rqmodel <- rq(value ~ cos_25C + sin_25C + decaycos + decaysin + decayint + cos_18C + sin_18C + 1,
                    data = temp_data, tau = 0.5)
      objective <- rqmodel$rho
      objective_storage[j] <- objective
      if (objective < objective_save) {
        ind_save <- j
        objective_save <- objective
        model_save <- rqmodel
      }
    }
  }
  guess <- as.vector(c(model_save$coefficients, lambda_range[ind_save]))
  return(guess)
}

####################
# function to generate weight for wild bootstrap, ref https://doi.org/10.1093/biomet/asr052
get_weight <- function(l = 73) {
  weight <- rep(NA, l)
  ind <- 0
  while (ind < l) {
    x <- runif(1, min = 3/4, max = 5/4)
    y <- runif(1, min = 0, max = 5/2)
    if (y <= 2*x) {
      ind <- ind + 1
      weight[ind] <- x
    }
  }
  rv <- sample(c(-1,1), size = l, replace = TRUE)
  weight <- weight * rv
  return(weight)
}

####################
# function to perform 1D&2D tests
test_2D <- function(Testdata) {
  test_mean <- as.matrix(apply(Testdata,2,mean))
  test_cov <- cov(Testdata)
  mDist <- rep(NA, nrow(Testdata))
  for (i in seq(nrow(Testdata))) {
    diff <- as.matrix(Testdata[i,]-test_mean)
    mDist[i] <- c(t(diff) %*% solve(test_cov) %*% diff)
  }
  mD <- c(t(test_mean) %*% solve(test_cov) %*% test_mean)
  test_pval <- sum(mDist>=mD)/nrow(Testdata)
  return(test_pval)
}

test_1D <- function(Testdata) {
  test_pval <- min(sum(Testdata<0),sum(Testdata>=0))*2/length(Testdata)
  return(test_pval)
}

####################
# function to perform wild bootstrap
wild_boot <- function(weight,fit_residual,fit_data,coarse,fine,boot_nretry,boot_fluc){
  data_boot <- data.frame(
    'value' = fit_data$value + fit_residual * weight,
    'time' = fit_data$time
  )
  par_est_boot <- get_guess(value = data_boot$value, time = data_boot$time,
                            coarse = seq(0.04,1,0.03), fine = seq(0.04,1,0.01))
  # get outputs
  ind <- 0
  errmsg <- TRUE
  while (errmsg & ind < boot_nretry) {
    ind <- ind + 1
    par_init <- par_est_boot * runif(9, min = 1-boot_fluc, max = 1+boot_fluc)
    bootmodel <- optim(par = par_init, fn = obj, gr = grad, method = 'L-BFGS-B',
                       value = data_boot$value, time = data_boot$time,
                       lower = c(rep(NA,8),0.03), upper = c(rep(NA,8),1.5))
    if (bootmodel$convergence == 0) {
      errmsg <- FALSE
    }
  }
  if (errmsg) {
    return(rep(NA,9))
  } else {
    return(bootmodel$par)
  }
}

####################
# function to get phase difference
get_angle <- function(a) {
  ang <- atan2(a[4],a[3]) - atan2(a[2],a[1])
  if (abs(ang) >= pi) {
    ang <- ang - 2*pi*sign(ang)
  }
  return(ang)
}

####################
# function to create plots
## full model plot
get_plots_nlrq <- function(result,data18C,data25C,time1,time2,time3) {
  TPM_18C <- as.vector(unlist(data18C))
  TPM_25C <- as.vector(unlist(data25C))
  temp_data <- data.frame(
    'TPM' = c(TPM_25C,TPM_18C),
    'time' = time1
  )
  temp_data_fitted <- data.frame(
    'time' = time2,
    'TPM' = func(c(unlist(result[3:11])),time2)
  )
  temp_data_extended <- data.frame(
    'time' = time3,
    'TPM' = func_type_positivetime(c(unlist(result[3:11])),time3,type=16,model='f')
  )
  anno1 <- sprintf('geneid = %s, genename = %s, lambda = %.4f',
                   result$geneid, result$gene_name, result$lambda)
  anno2 <- sprintf('pval_cyc25C = %.3f, pval_cyc18C = %.3f, pval_cycdiff = %.3f',
                   result$pval_cycle_25C, result$pval_cycle_18C, result$pval_cycle_diff)
  anno3 <- sprintf('pval_phasediff = %.3f, pval_ampdiff = %.3f, pval_mediandiff = %.3f',
                   result$pval_phase_diff, result$pval_amp_diff, result$pval_median_diff)
  anno4 <- sprintf('amp25C = %.2f, phase25C = %.2f, amp18C = %.2f, phase18C = %.2f',
                   sqrt(result$sin_25C^2+result$cos_25C^2), 
                   atan2(result$sin_25C,result$cos_25C),
                   sqrt((result$sin_18C+result$sin_25C)^2+(result$cos_18C+result$cos_25C)^2), 
                   atan2(result$sin_18C+result$sin_25C,result$cos_18C+result$cos_25C))
  g <- ggplot(data = temp_data[1:sum(time1<0),], aes(time,TPM)) + geom_point(size = 4, color = 'red') +
    geom_point(data = temp_data[(sum(time1<0)+1):length(time1),], aes(time,TPM), size = 4, color = 'blue') +
    geom_line(data = temp_data_fitted, aes(time,TPM), color = 'black', lwd = 1.5) +
    geom_line(data = temp_data_extended, aes(time,TPM), color = 'red', lwd = 1, linetype = "twodash") +
    ylim(min(c(temp_data$TPM,temp_data_extended$TPM,temp_data_fitted$TPM))-5,NA) +
    scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
    geom_text(aes(x=-Inf,y=Inf,hjust=-0.1,vjust=3,label=anno1), size = 3) +
    geom_text(aes(x=-Inf,y=Inf,hjust=-0.1,vjust=5,label=anno2), size = 3) +
    geom_text(aes(x=-Inf,y=Inf,hjust=-0.1,vjust=7,label=anno3), size = 3) +
    geom_text(aes(x=-Inf,y=Inf,hjust=-0.1,vjust=9,label=anno4), size = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
          panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
  return(g)
}

# plot for other types
get_plots_nlrq_type <- function(par,data18C,data25C,time1,time2,time3,model,type,lambda) {
  modelnum <- as.numeric(substr(model,1,2))
  modelchar <- substr(model,3,3)
  TPM_18C <- as.vector(unlist(data18C))
  TPM_25C <- as.vector(unlist(data25C))
  temp_data <- data.frame(
    'TPM' = c(TPM_25C,TPM_18C),
    'time' = time1
  )
  if (modelchar == "f") {
    temp_data_fitted <- data.frame(
      'time' = time2,
      'TPM' = func_fullmodel_type(par=par,time=time2,type=modelnum)
    )
  } else if (modelchar == "p") {
    temp_data_fitted <- data.frame(
      'time' = time2,
      'TPM' = func_partialmodel_type(par=par,time=time2,type=modelnum)
    )
  }
  temp_data_extended <- data.frame(
    'time' = time3,
    'TPM' = func_type_positivetime(par=par,time=time3,type=modelnum,model=modelchar)
  )
  
  anno1 <- sprintf('model = %s, type = %i, lambda = %.4f', model, type, lambda)
  g <- ggplot(data = temp_data[1:sum(time1<0),], aes(time,TPM)) + geom_point(size = 4, color = 'red') +
    geom_point(data = temp_data[(sum(time1<0)+1):length(time1),], aes(time,TPM), size = 4, color = 'blue') +
    geom_line(data = temp_data_fitted, aes(time,TPM), color = 'black', lwd = 1.5) +
    geom_line(data = temp_data_extended, aes(time,TPM), color = 'red', lwd = 1, linetype = "twodash") +
    ylim(min(c(temp_data$TPM,temp_data_extended$TPM,temp_data_fitted$TPM))-5,NA) +
    scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
    geom_text(aes(x=-Inf,y=Inf,hjust=-0.1,vjust=3,label=anno1), size = 3) +
    theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
          panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
  return(g)
}

####################
# function to create plots, no labels
## full model plot
get_plots_nlrq_nolabel <- function(result,data18C,data25C,time1,time2,time3) {
  TPM_18C <- as.vector(unlist(data18C))
  TPM_25C <- as.vector(unlist(data25C))
  temp_data <- data.frame(
    'TPM' = c(TPM_25C,TPM_18C),
    'time' = time1
  )
  temp_data_fitted <- data.frame(
    'time' = time2,
    'TPM' = func(c(unlist(result[3:11])),time2)
  )
  temp_data_extended <- data.frame(
    'time' = time3,
    'TPM' = func_type_positivetime(c(unlist(result[3:11])),time3,type=16,model='f')
  )
  g <- ggplot(data = temp_data[1:sum(time1<0),], aes(time,TPM)) + geom_point(size = 4, color = 'red') +
    geom_point(data = temp_data[(sum(time1<0)+1):length(time1),], aes(time,TPM), size = 4, color = 'blue') +
    geom_line(data = temp_data_fitted, aes(time,TPM), color = 'black', lwd = 1.5) +
    geom_line(data = temp_data_extended, aes(time,TPM), color = 'red', lwd = 1, linetype = "twodash") +
    ylim(min(c(temp_data$TPM,temp_data_extended$TPM,temp_data_fitted$TPM))-5,NA) +
    scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
    theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
          panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
  return(g)
}

# plot for other types
get_plots_nlrq_type_nolabel <- function(par,data18C,data25C,time1,time2,time3,model,type,lambda) {
  modelnum <- as.numeric(substr(model,1,2))
  modelchar <- substr(model,3,3)
  TPM_18C <- as.vector(unlist(data18C))
  TPM_25C <- as.vector(unlist(data25C))
  temp_data <- data.frame(
    'TPM' = c(TPM_25C,TPM_18C),
    'time' = time1
  )
  if (modelchar == "f") {
    temp_data_fitted <- data.frame(
      'time' = time2,
      'TPM' = func_fullmodel_type(par=par,time=time2,type=modelnum)
    )
  } else if (modelchar == "p") {
    temp_data_fitted <- data.frame(
      'time' = time2,
      'TPM' = func_partialmodel_type(par=par,time=time2,type=modelnum)
    )
  }
  temp_data_extended <- data.frame(
    'time' = time3,
    'TPM' = func_type_positivetime(par=par,time=time3,type=modelnum,model=modelchar)
  )
  
  g <- ggplot(data = temp_data[1:sum(time1<0),], aes(time,TPM)) + geom_point(size = 4, color = 'red') +
    geom_point(data = temp_data[(sum(time1<0)+1):length(time1),], aes(time,TPM), size = 4, color = 'blue') +
    geom_line(data = temp_data_fitted, aes(time,TPM), color = 'black', lwd = 1.5) +
    geom_line(data = temp_data_extended, aes(time,TPM), color = 'red', lwd = 1, linetype = "twodash") +
    ylim(min(c(temp_data$TPM,temp_data_extended$TPM,temp_data_fitted$TPM))-5,NA) +
    scale_x_continuous(breaks = seq(-120, 120, by = 24)) +
    theme(panel.background = element_rect(fill = "white", colour = "black", linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"),
          panel.grid.minor = element_line(size = 0.25, linetype = "solid", colour = "grey"))
  return(g)
}

####################
# function to compute PRC
get_PRC <- function(A1,A2,B1,B2) {
  phi1 <- atan2(A2,A1)
  phi2 <- atan2(A2+B2,A1+B1)
  if (phi1 > 0) {
    leftcomp <- phi1/pi*12
  } else {
    leftcomp <- (2*pi+phi1)/pi*12
  }
  if (phi2 > 0) {
    rightcomp <- phi2/pi*12
  } else {
    rightcomp <- (2*pi+phi2)/pi*12
  }
  T0 <- 24
  T1 <- leftcomp + rightcomp
  result <- c(leftcomp,rightcomp,T1/T0-1)
  return(result)
}

####################
# function for all types of full model, used in optim
func_fullmodel_type <- function(par,time,type) {
  result <- time * 0
  if (type == 1) {
    
    result <- par[1] +
      par[2] * exp(-par[4]*time)*(cos(time*pi/12)-1)*(time>=0) +
      par[3] * exp(-par[4]*time)*sin(time*pi/12)*(time>=0)
    
  } else if (type == 2) {
    
    result <- par[1] +
      par[2] * (exp(-par[5]*time)*cos(time*pi/12)-1)*(time>=0) +
      par[3] * exp(-par[5]*time)*sin(time*pi/12)*(time>=0) +
      par[4] * (exp(-par[5]*time)-1)*(time>=0)
    
  } else if (type == 7) {
    
    result <- par[1] +
      par[2] * (exp(-par[6]*time)-1)*cos(time*pi/12)*(time>=0) +
      par[3] * exp(-par[6]*time)*sin(time*pi/12)*(time>=0) +
      par[4] * (exp(-par[6]*time)-cos(time*pi/12))*(time>=0) +
      par[5] * sin(time*pi/12)*(time>=0)
    
  } else if (type == 8) {
    
    result <- par[1] +
      par[2] * (exp(-par[7]*time)*cos(time*pi/12)-1)*(time>=0) +
      par[3] * exp(-par[7]*time)*sin(time*pi/12)*(time>=0) +
      par[4] * (exp(-par[7]*time)-1)*(time>=0) +
      par[5] * (cos(time*pi/12)-1)*(time>=0) +
      par[6] * sin(time*pi/12)*(time>=0)
    
  } else if (type == 11) {
    
    result <- par[1] + 
      par[2] * (cos(time*pi/12)*(time<0)+exp(-par[6]*time)*(time>=0)) + 
      par[3] * sin(time*pi/12)*(time<0) +
      par[4] * exp(-par[6]*time)*(cos(time*pi/12)-1)*(time>=0) +
      par[5] * exp(-par[6]*time)*sin(time*pi/12)*(time>=0)
    
  } else if (type == 12) {
    
    result <- par[1] + 
      par[2] * (cos(time*pi/12)*(time<0)+(time>=0)) +
      par[3] * sin(time*pi/12)*(time<0) +
      par[4] * (exp(-par[7]*time)*cos(time*pi/12)-1)*(time>=0) +
      par[5] * exp(-par[7]*time)*sin(time*pi/12)*(time>=0) +
      par[6] * (exp(-par[7]*time)-1)*(time>=0)
    
  } else if (type == 13) {
    
    result <- par[1] + 
      par[2] * cos(time*pi/12) + 
      par[3] * sin(time*pi/12) +
      par[4] * exp(-par[6]*time)*(cos(time*pi/12)-1)*(time>=0) +
      par[5] * exp(-par[6]*time)*sin(time*pi/12)*(time>=0)
    
  } else if (type == 14) {
    
    result <- par[1] + 
      par[2] * cos(time*pi/12) + 
      par[3] * sin(time*pi/12) +
      par[4] * (exp(-par[7]*time)*cos(time*pi/12)-1)*(time>=0) +
      par[5] * exp(-par[7]*time)*sin(time*pi/12)*(time>=0) +
      par[6] * (exp(-par[7]*time)-1)*(time>=0)
    
  } else if (type == 15) {
    
    result <- par[1] + 
      par[2] * cos(time*pi/12) + 
      par[3] * sin(time*pi/12) +
      par[4] * (exp(-par[8]*time)-1)*cos(time*pi/12)*(time>=0) +
      par[5] * exp(-par[8]*time)*sin(time*pi/12)*(time>=0) +
      par[6] * (exp(-par[8]*time)-cos(time*pi/12))*(time>=0) +
      par[7] * sin(time*pi/12)*(time>=0)
    
  } else if (type == 16) {
    
    result <- func(par, time)
    
  }
  return(result)
}

####################
# objective function to minimize all types of full model, used in optim
obj_fullmodel_type <- function(par,value,time,type) {
  AD <- value - func_fullmodel_type(par,time,type)
  result <- 1/2 * sum(AD*sign(AD))
  return(result)
}

####################
# gradient of all types of full model's objective function, used in optim
grad_fullmodel_type <- function(par,value,time,type) {
  abs_fac <- sign(value - func_fullmodel_type(par,time,type))
  if (type == 1) {
    
    result <- rep(NA,4)
    result[1] <- -1/2 * sum(abs_fac * 1)
    result[2] <- -1/2 * sum(abs_fac * exp(-par[4]*time)*(cos(time*pi/12)-1)*(time>=0))
    result[3] <- -1/2 * sum(abs_fac * exp(-par[4]*time)*sin(time*pi/12)*(time>=0))
    result[4] <- -1/2 * sum(abs_fac * -time*exp(-par[4]*time)*(par[2]*cos(time*pi/12)+par[3]*sin(time*pi/12)-par[2])
                            *(time>=0))
    
  } else if (type == 2) {
    
    result <- rep(NA,5)
    result[1] <- -1/2 * sum(abs_fac * 1)
    result[2] <- -1/2 * sum(abs_fac * (exp(-par[5]*time)*cos(time*pi/12)-1)*(time>=0))
    result[3] <- -1/2 * sum(abs_fac * (exp(-par[5]*time)*sin(time*pi/12)-0)*(time>=0))
    result[4] <- -1/2 * sum(abs_fac * (exp(-par[5]*time)-1)*(time>=0))
    result[5] <- -1/2 * sum(abs_fac * -time*exp(-par[5]*time)*(par[2]*cos(time*pi/12)+par[3]*sin(time*pi/12)+par[4])
                            *(time>=0))
    
  } else if (type == 7) {
    
    result <- rep(NA,6)
    result[1] <- -1/2 * sum(abs_fac * 1)
    result[2] <- -1/2 * sum(abs_fac * (exp(-par[6]*time)-1)*cos(time*pi/12)*(time>=0))
    result[3] <- -1/2 * sum(abs_fac * exp(-par[6]*time)*sin(time*pi/12)*(time>=0))
    result[4] <- -1/2 * sum(abs_fac * (exp(-par[6]*time)-cos(time*pi/12))*(time>=0))
    result[5] <- -1/2 * sum(abs_fac * sin(time*pi/12)*(time>=0))
    result[6] <- -1/2 * sum(abs_fac * -time*exp(-par[6]*time)*
                              (par[2]*cos(time*pi/12)+par[3]*sin(time*pi/12)+par[4])*(time>=0))
    
  } else if (type == 8) {
    
    result <- rep(NA,7)
    result[1] <- -1/2 * sum(abs_fac * 1)
    result[2] <- -1/2 * sum(abs_fac * (exp(-par[7]*time)*cos(time*pi/12)-1)*(time>=0))
    result[3] <- -1/2 * sum(abs_fac * (exp(-par[7]*time)*sin(time*pi/12)-0)*(time>=0))
    result[4] <- -1/2 * sum(abs_fac * (exp(-par[7]*time)-1)*(time>=0))
    result[5] <- -1/2 * sum(abs_fac * (cos(time*pi/12)-1)*(time>=0))
    result[6] <- -1/2 * sum(abs_fac * sin(time*pi/12)*(time>=0))
    result[7] <- -1/2 * sum(abs_fac * -time*exp(-par[7]*time)*(par[2]*cos(time*pi/12)+par[3]*sin(time*pi/12)+par[4])
                            *(time>=0))
    
  } else if (type == 11) {
    
    result <- rep(NA,6)
    result[1] <- -1/2 * sum(abs_fac * 1)
    result[2] <- -1/2 * sum(abs_fac * (cos(time*pi/12)*(time<0) + exp(-par[6]*time)*(time>=0)))
    result[3] <- -1/2 * sum(abs_fac * sin(time*pi/12)*(time<0))
    result[4] <- -1/2 * sum(abs_fac * exp(-par[6]*time)*(cos(time*pi/12)-1)*(time>=0))
    result[5] <- -1/2 * sum(abs_fac * exp(-par[6]*time)*sin(time*pi/12)*(time>=0))
    result[6] <- -1/2 * sum(abs_fac * -time*exp(-par[6]*time)*
                              (par[4]*(cos(time*pi/12)-1)+par[5]*sin(time*pi/12)+par[2])*(time>=0))
    
  } else if (type == 12) {
    
    result <- rep(NA,7)
    result[1] <- -1/2 * sum(abs_fac * 1)
    result[2] <- -1/2 * sum(abs_fac * (cos(time*pi/12)*(time<0)+(time>=0)))
    result[3] <- -1/2 * sum(abs_fac * sin(time*pi/12)*(time<0))
    result[4] <- -1/2 * sum(abs_fac * (exp(-par[7]*time)*cos(time*pi/12)-1)*(time>=0))
    result[5] <- -1/2 * sum(abs_fac * (exp(-par[7]*time)*sin(time*pi/12)-0)*(time>=0))
    result[6] <- -1/2 * sum(abs_fac * (exp(-par[7]*time)-1)*(time>=0))
    result[7] <- -1/2 * sum(abs_fac * -time*exp(-par[7]*time)*(par[4]*cos(time*pi/12)+par[5]*sin(time*pi/12)+par[6])
                            *(time>=0))
    
  } else if (type == 13) {
    
    result <- rep(NA,6)
    result[1] <- -1/2 * sum(abs_fac * 1)
    result[2] <- -1/2 * sum(abs_fac * cos(time*pi/12))
    result[3] <- -1/2 * sum(abs_fac * sin(time*pi/12))
    result[4] <- -1/2 * sum(abs_fac * exp(-par[6]*time)*(cos(time*pi/12)-1)*(time>=0))
    result[5] <- -1/2 * sum(abs_fac * exp(-par[6]*time)*sin(time*pi/12)*(time>=0))
    result[6] <- -1/2 * sum(abs_fac * -time*exp(-par[6]*time)*(par[4]*cos(time*pi/12)+par[5]*sin(time*pi/12)-par[4])
                            *(time>=0))
    
  } else if (type == 14) {
    
    result <- rep(NA,7)
    result[1] <- -1/2 * sum(abs_fac * 1)
    result[2] <- -1/2 * sum(abs_fac * cos(time*pi/12))
    result[3] <- -1/2 * sum(abs_fac * sin(time*pi/12))
    result[4] <- -1/2 * sum(abs_fac * (exp(-par[7]*time)*cos(time*pi/12)-1)*(time>=0))
    result[5] <- -1/2 * sum(abs_fac * (exp(-par[7]*time)*sin(time*pi/12)-0)*(time>=0))
    result[6] <- -1/2 * sum(abs_fac * (exp(-par[7]*time)-1)*(time>=0))
    result[7] <- -1/2 * sum(abs_fac * -time*exp(-par[7]*time)*(par[4]*cos(time*pi/12)+par[5]*sin(time*pi/12)+par[6])
                            *(time>=0))
    
  } else if (type == 15) {
    
    result <- rep(NA,8)
    result[1] <- -1/2 * sum(abs_fac * 1)
    result[2] <- -1/2 * sum(abs_fac * cos(time*pi/12))
    result[3] <- -1/2 * sum(abs_fac * sin(time*pi/12))
    result[4] <- -1/2 * sum(abs_fac * (exp(-par[8]*time)-1)*cos(time*pi/12)*(time>=0))
    result[5] <- -1/2 * sum(abs_fac * exp(-par[8]*time)*sin(time*pi/12)*(time>=0))
    result[6] <- -1/2 * sum(abs_fac * (exp(-par[8]*time)-cos(time*pi/12))*(time>=0))
    result[7] <- -1/2 * sum(abs_fac * sin(time*pi/12)*(time>=0))
    result[8] <- -1/2 * sum(abs_fac * -time*exp(-par[8]*time)*
                              (par[4]*cos(time*pi/12)+par[5]*sin(time*pi/12)+par[6])*(time>=0))
    
  }
  return(result)
}

####################
# function to generate an initial parameter guess using sweep, for all types of model
get_guess_fullmodel_type <- function(value,time,coarse=seq(0.04,1,0.03),fine=NULL,type) {
  lambda_range <- coarse
  objective_storage <- 0*lambda_range
  objective_save <- Inf
  ind_save <- 0
  model_save <- NA
  for (j in seq(length(lambda_range))) {
    lambda <- lambda_range[j]
    temp_data <- data.frame(
      'value' = value,
      'time' = time,
      'cos_25C' = cos(time*pi/12),
      'sin_25C' = sin(time*pi/12),
      'decaycos' = (exp(-lambda*time)*cos(time*pi/12)-1)*(time>=0),
      'decaysin' = exp(-lambda*time)*sin(time*pi/12)*(time>=0),
      'decayint' = (exp(-lambda*time)-1)*(time>=0),
      'decaycos_add' = exp(-lambda*time)*(cos(time*pi/12)-1)*(time>=0),
      'decaycos_add2' = (exp(-lambda*time)-1)*cos(time*pi/12)*(time>=0),
      'decayint_add' = (exp(-lambda*time)-cos(time*pi/12))*(time>=0),
      'cos_18C' = (cos(time*pi/12)-1)*(time>=0),
      'sin_18C' = sin(time*pi/12)*(time>=0),
      'cos_25C_add' = cos(time*pi/12)*(time<0)+(time>=0),
      'cos_25C_add2' = (cos(time*pi/12)*(time<0)+exp(-lambda*time)*(time>=0)),
      'sin_25C_add' = sin(time*pi/12)*(time<0)
    )
    if (type == 1) {
      rqmodel <- rq(value ~ decaycos_add + decaysin + 1,
                    data = temp_data, tau = 0.5)
    } else if (type == 2) {
      rqmodel <- rq(value ~ decaycos + decaysin + decayint + 1,
                    data = temp_data, tau = 0.5)
    } else if (type == 7) {
      rqmodel <- rq(value ~ decaycos_add2 + decaysin + decayint_add + sin_18C + 1,
                    data = temp_data, tau = 0.5)
    } else if (type == 8) {
      rqmodel <- rq(value ~ decaycos + decaysin + decayint + cos_18C + sin_18C + 1,
                    data = temp_data, tau = 0.5)
    } else if (type == 11) {
      rqmodel <- rq(value ~ cos_25C_add2 + sin_25C_add + decaycos_add + decaysin + 1,
                    data = temp_data, tau = 0.5)
    } else if (type == 12) {
      rqmodel <- rq(value ~ cos_25C_add + sin_25C_add + decaycos + decaysin + decayint + 1,
                    data = temp_data, tau = 0.5)
    } else if (type == 13) {
      rqmodel <- rq(value ~ cos_25C + sin_25C + decaycos_add + decaysin + 1,
                    data = temp_data, tau = 0.5)
    } else if (type == 14) {
      rqmodel <- rq(value ~ cos_25C + sin_25C + decaycos + decaysin + decayint + 1,
                    data = temp_data, tau = 0.5)
    } else if (type == 15) {
      rqmodel <- rq(value ~ cos_25C + sin_25C + decaycos_add2 + decaysin + decayint_add + sin_18C + 1,
                    data = temp_data, tau = 0.5)
    }
    objective <- rqmodel$rho
    objective_storage[j] <- objective
    if (objective < objective_save) {
      ind_save <- j
      objective_save <- objective
      model_save <- rqmodel
    }
  }
  if (!is.null(fine) & (lambda_range[ind_save] == max(lambda_range) || lambda_range[ind_save] == min(lambda_range))) {
    lambda_range <- fine
    objective_storage <- 0*lambda_range
    objective_save <- Inf
    ind_save <- 0
    model_save <- NA
    for (j in seq(length(lambda_range))) {
      lambda <- lambda_range[j]
      temp_data <- data.frame(
        'value' = value,
        'time' = time,
        'cos_25C' = cos(time*pi/12),
        'sin_25C' = sin(time*pi/12),
        'decaycos' = (exp(-lambda*time)*cos(time*pi/12)-1)*(time>=0),
        'decaysin' = exp(-lambda*time)*sin(time*pi/12)*(time>=0),
        'decayint' = (exp(-lambda*time)-1)*(time>=0),
        'decaycos_add' = exp(-lambda*time)*(cos(time*pi/12)-1)*(time>=0),
        'decaycos_add2' = (exp(-lambda*time)-1)*cos(time*pi/12)*(time>=0),
        'decayint_add' = (exp(-lambda*time)-cos(time*pi/12))*(time>=0),
        'cos_18C' = (cos(time*pi/12)-1)*(time>=0),
        'sin_18C' = sin(time*pi/12)*(time>=0),
        'cos_25C_add' = cos(time*pi/12)*(time<0)+(time>=0),
        'cos_25C_add2' = (cos(time*pi/12)*(time<0)+exp(-lambda*time)*(time>=0)),
        'sin_25C_add' = sin(time*pi/12)*(time<0)
      )
      if (type == 1) {
        rqmodel <- rq(value ~ decaycos_add + decaysin + 1,
                      data = temp_data, tau = 0.5)
      } else if (type == 2) {
        rqmodel <- rq(value ~ decaycos + decaysin + decayint + 1,
                      data = temp_data, tau = 0.5)
      } else if (type == 7) {
        rqmodel <- rq(value ~ decaycos_add2 + decaysin + decayint_add + sin_18C + 1,
                      data = temp_data, tau = 0.5)
      } else if (type == 8) {
        rqmodel <- rq(value ~ decaycos + decaysin + decayint + cos_18C + sin_18C + 1,
                      data = temp_data, tau = 0.5)
      } else if (type == 11) {
        rqmodel <- rq(value ~ cos_25C_add2 + sin_25C_add + decaycos_add + decaysin + 1,
                      data = temp_data, tau = 0.5)
      } else if (type == 12) {
        rqmodel <- rq(value ~ cos_25C_add + sin_25C_add + decaycos + decaysin + decayint + 1,
                      data = temp_data, tau = 0.5)
      } else if (type == 13) {
        rqmodel <- rq(value ~ cos_25C + sin_25C + decaycos_add + decaysin + 1,
                      data = temp_data, tau = 0.5)
      } else if (type == 14) {
        rqmodel <- rq(value ~ cos_25C + sin_25C + decaycos + decaysin + decayint + 1,
                      data = temp_data, tau = 0.5)
      } else if (type == 15) {
        rqmodel <- rq(value ~ cos_25C + sin_25C + decaycos_add2 + decaysin + decayint_add + sin_18C + 1,
                      data = temp_data, tau = 0.5)
      }
      objective <- rqmodel$rho
      objective_storage[j] <- objective
      if (objective < objective_save) {
        ind_save <- j
        objective_save <- objective
        model_save <- rqmodel
      }
    }
  }
  guess <- as.vector(c(model_save$coefficients, lambda_range[ind_save]))
  return(guess)
}

####################
# functions fo get AICc
## partial models of all types
get_partialmodel_AICc <- function(data,type) {
  n <- length(data$time)
  temp_data <- data.frame(
    'value' = data$value,
    'time' = data$time,
    'cos_25C' = cos(data$time*pi/12),
    'sin_25C' = sin(data$time*pi/12),
    'cos_18C' = (cos(data$time*pi/12)-1)*(data$time>=0),
    'sin_18C' = (sin(data$time*pi/12)-0)*(data$time>=0),
    'ind_time' = as.numeric(data$time>=0),
    'cos_25C_add' = cos(data$time*pi/12)*(data$time<0)+(data$time>=0),
    'sin_25C_add' = sin(data$time*pi/12)*(data$time<0)
  )
  if (type == 1) {
    partialmodel <- rq(value ~ 1, data = temp_data)
    k <- 1
    result2 <- 0
  } else if (type == 2) {
    partialmodel <- rq(value ~ ind_time + 1, data = temp_data)
    k <- 2
    result2 <- Inf
  } else if (type == 7) {
    partialmodel <- rq(value ~ sin_18C + 1, data = temp_data)
    k <- 2
    result2 <- Inf
  } else if (type == 8) {
    partialmodel <- rq(value ~ cos_18C + sin_18C + 1, data = temp_data)
    k <- 3
    result2 <- Inf
  } else if (type == 11) {
    partialmodel <- rq(value ~ sin_25C_add + 1, data = temp_data)
    k <- 2
    result2 <- Inf
  } else if (type == 12) {
    partialmodel <- rq(value ~ cos_25C_add + sin_25C_add + 1, data = temp_data)
    k <- 3
    result2 <- Inf
  } else if (type == 13) {
    partialmodel <- rq(value ~ cos_25C + sin_25C + 1, data = temp_data)
    k <- 3
    result2 <- 0
  } else if (type == 14) {
    partialmodel <- rq(value ~ cos_25C + sin_25C + ind_time + 1, data = temp_data)
    k <- 4
    result2 <- Inf
  } else if (type == 15) {
    partialmodel <- rq(value ~ cos_25C + sin_25C + sin_18C + 1, data = temp_data)
    k <- 4
    result2 <- Inf
  } else if (type == 16) {
    partialmodel <- rq(value ~ cos_25C + sin_25C + cos_18C + sin_18C + 1, data = temp_data)
    k <- 5
    result2 <- Inf
  }
  result1 <- 2*k - 2*(n*log(n/4/partialmodel$rho)-n) + (2*k^2+2*k)/(n-k-1)
  result3 <- partialmodel$coefficients
  return(list(result1,result2,result3))
}

####################
## full models of all types
get_fullmodel_AICc <- function(data,nretry,type) {
  n <- length(data$time)
  par_est_type <- get_guess_fullmodel_type(value = data$value, time = data$time,
                                           coarse = seq(0.04,1,0.03),
                                           fine = seq(0.04,1,0.01), type = type)
  if (type == 1) {
    k <- 4
  } else if (type == 2) {
    k <- 5
  } else if (type == 7) {
    k <- 6
  } else if (type == 8) {
    k <- 7
  } else if (type == 11) {
    k <- 6
  } else if (type == 12) {
    k <- 7
  } else if (type == 13) {
    k <- 6
  } else if (type == 14) {
    k <- 7
  } else if (type == 15) {
    k <- 8
  }
  ind <- 0
  errmsg <- TRUE
  while (errmsg & ind < nretry) {
    ind <- ind + 1
    par_init_type <- par_est_type * runif(k, min = 1-fluc, max = 1+fluc)
    nlrqmodel <- optim(par = par_init_type, fn = obj_fullmodel_type, 
                       gr = grad_fullmodel_type, method = 'L-BFGS-B',
                       value = data$value, time = data$time, type = type,
                       lower = c(rep(NA,k-1),0.03), upper = c(rep(NA,k-1),1.5))
    if (nlrqmodel$convergence == 0) {
      errmsg <- FALSE
    }
  }
  result1 <- 2*k - 2*(n*log(n/4/nlrqmodel$value)-n) + (2*k^2+2*k)/(n-k-1)
  result2 <- nlrqmodel$par[k]
  result3 <- nlrqmodel$par
  return(list(result1,result2,result3))
}

####################
# function for all types of partial model
func_partialmodel_type <- function(par,time,type) {
  result <- time * 0
  if (type == 1) {
    
    result <- rep(par[1],length(time))
    
  } else if (type == 2) {
    
    result <- par[1] + 
      par[2] * (time>=0)
    
  } else if (type == 7) {
    
    result <- par[1] + 
      par[2] * sin(time*pi/12) * (time>=0)
    
  } else if (type == 8) {
    
    result <- par[1] + 
      par[2] * (cos(time*pi/12)-1) * (time>=0) +
      par[3] * sin(time*pi/12) * (time>=0)
    
  } else if (type == 11) {
    
    result <- par[1] + 
      par[2] * sin(time*pi/12) * (time<0)
    
  } else if (type == 12) {
    
    result <- par[1] + 
      par[2] * (cos(time*pi/12) * (time<0) + (time>=0))+ 
      par[3] * sin(time*pi/12) * (time<0)
    
  } else if (type == 13) {
    
    result <- par[1] + 
      par[2] * cos(time*pi/12) + 
      par[3] * sin(time*pi/12)
    
  } else if (type == 14) {
    
    result <- par[1] + 
      par[2] * cos(time*pi/12) + 
      par[3] * sin(time*pi/12) +
      par[4] * (time>=0)
    
  } else if (type == 15) {
    
    result <- par[1] + 
      par[2] * cos(time*pi/12) + 
      par[3] * sin(time*pi/12) +
      par[4] * sin(time*pi/12) * (time>=0)
    
  } else if (type == 16) {
    
    result <- par[1] + 
      par[2] * cos(time*pi/12) +
      par[3] * sin(time*pi/12) +
      par[4] * (cos(time*pi/12)-1) * (time>=0) +
      par[5] * sin(time*pi/12) * (time>=0)
    
  }
  return(result)
}

# Get subset of fitted results for a gene id list (in a GO category)
results_from_FlybaseIDs <- function(IDList,result1,result2) {
  result <- matrix(NA,nrow=length(IDList),ncol=16)
  colnames(result) <- c('geneid','genename','geneindex','pval_cycle_25C','pval_cycle_18C','pval_cycle_diff',
                        'pval_phase_diff','pval_amp_diff','pval_median_diff','type','model',
                        'TRUE_lambda','typematch','TransSig','dAICc','nextmodel')
  result <- as.data.frame(result)
  colvec <- c('01f','02f','07f','08f','11f','12f','13f','14f','15f','16f',
              '01p','02p','07p','08p','11p','12p','13p','14p','15p','16p')
  for (i in seq(IDList)) {
    j <- which(result1$geneid == IDList[i])
    if (length(j) > 0) {
      result[i,1:2] <- result1[j,1:2]
      result[i,4:9] <- result1[j,12:17]
      result$type[i] <- result2$type[j]
      result$model[i] <- result2$model[j]
      result$TRUE_lambda[i] <- result2$TRUE_lambda[j]
      result$geneindex[i] <- j
      AICc_vec <- as.vector(unlist(result2[j,4:23]))
      result$nextmodel[i] <- colvec[match(sort(AICc_vec)[2],AICc_vec)]
      AICc_vec <- sort(AICc_vec)
      result$dAICc[i] <- diff(AICc_vec)[1]
      if (!is.na(result$model[i])) {
        if (result$type[i] == as.numeric(substr(result$model[i],1,2))) {
          result$typematch[i] <- 1
        } else {
          result$typematch[i] <- 0
        }
        if (substr(result$model[i],3,3) == "f") {
          result$TransSig[i] <- 1
        } else {
          result$TransSig[i] <- 0
        }
      }
    }
  }
  result <- result[order(result$typematch,result$TransSig,result$TRUE_lambda,
                         result$model,decreasing=TRUE),]
  return(result)
}

####################
# function to produce curves used in plotting, with the curve at positive be an extension of negative time
func_type_positivetime <- function(par,time,type,model) {
  result <- time * 0
  if (type == 1) {
    result <- rep(par[1],length(time))
  } else if (type == 2) {
    result <- rep(par[1],length(time))
  } else if (type == 7) {
    result <- rep(par[1],length(time))
  } else if (type == 8) {
    result <- rep(par[1],length(time))
  } else if (type == 11 & model == 'p') {
    result <- par[1] + par[2] * sin(time*pi/12)
  } else if (type == 11 & model == 'f') {
    result <- par[1] + par[2] * cos(time*pi/12) + par[3] * sin(time*pi/12)
  } else if (type == 12) {
    result <- par[1] + par[2] * cos(time*pi/12) + par[3] * sin(time*pi/12)
  } else if (type == 13) {
    result <- par[1] + par[2] * cos(time*pi/12) + par[3] * sin(time*pi/12)
  } else if (type == 14) {
    result <- par[1] + par[2] * cos(time*pi/12) + par[3] * sin(time*pi/12)
  } else if (type == 15) {
    result <- par[1] + par[2] * cos(time*pi/12) + par[3] * sin(time*pi/12)
  } else if (type == 16) {
    result <- par[1] + par[2] * cos(time*pi/12) + par[3] * sin(time*pi/12)
  }
  return(result)
}