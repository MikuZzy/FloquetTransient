# This is the source function file specifically for network analysis. 

# Function to get lambdas from gene ids.
lambda_from_geneid <- function(geneid,results) {
  ind <- which(results$geneid==geneid)
  return(results$TRUE_lambda[ind])
}

# Function to get gene names from gene ids.
genename_from_geneid <- function(geneid,results) {
  ind <- which(results$geneid==geneid)
  return(results$gene_name[ind])
}

# Function to get gene ids from gene names.
geneid_from_genename <- function(genename,results) {
  ind <- which(results$gene_name==genename)
  return(results$geneid[ind])
}

# Function to get lambdas from gene names.
lambda_from_genename <- function(genename,results) {
  ind <- which(results$gene_name==genename)
  return(results$TRUE_lambda[ind])
}

# Function to get best models from gene names.
model_from_genename <- function(genename,results) {
  ind <- which(results$gene_name==genename)
  return(results$model[ind])
}

# Function to get detected cycling status from gene names.
cycling_from_genename <- function(genename,results,temp) {
  if (temp==18) {
    vec <- c("07","08","13","14","15","16")
  } else if (temp==25) {
    vec <- c("11","12","13","14","15","16")
  }
  ind <- which(results$gene_name==genename)
  model <- results$model[ind]
  if (substr(model,1,2) %in% vec) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Function to generate directed graph from a data frame.
get_graph <- function(graph_df) {
  df <- data.frame(
    'Node1' <- graph_df$src,
    'Node2' <- graph_df$dest
  )
  Dir_Graph <- graph.data.frame(df,directed = TRUE, vertices = NULL)
  return(Dir_Graph)
}

# Function to generate customized graph layout based on fitted lambdas. Genes are divided into different 
# chunks based on percentiles of their lambdas. The vertical axis correspond to ranks of lambda, and the 
# horizontal axis within each chunk are evenly distributed and optimized (by randomly switch pairs of 
# horizontal coordinates within each chunks) so that the total plot length of edges are minimized.
#
# Input parameters:
#   graph: igraph object to be plotted
#   lambda: vector of fitted lambdas
#   xlim: 2-element vector that specify plot ranges for horizontal axis
#   ylim: 2-element vector that specify plot ranges for vertical axis
#   breaks: vector of percentiles to divide genes into different chunks
#   ntrials: number of trials to perform randomly switching of paris of horizontal coordinates
#   counts: a normalization parameter, together with `tol`, to scale the xlim within each chunk
#   tol: a normalization parameter, together with `counts`, to scale the xlim within each chunk
#
# Output:
#   result: 2-column matrix of coordinates as the plot layout
get_layout <- function(graph,lambda,xlim,ylim,breaks,ntrial,counts,tol) {
  n <- length(breaks)
  result <- matrix(0,nrow=gorder(graph),ncol=2)
  result[,2] <- (rank(lambda,ties.method="max")-1)/(gorder(graph)-1)*(ylim[2]-ylim[1])+ylim[1]
  grp_info <- vector("list",n)
  
  # assign initial horizontal coordinates within each chunk
  for (i in seq(n-1)) {
    grp_info[[i]] <- which(lambda>=breaks[i] & lambda<breaks[i+1])
    if (length(grp_info[[i]])/counts>tol) {
      xgrid <- seq(0,xlim[2]-xlim[1],length.out=length(grp_info[[i]])+2)*length(grp_info[[i]])/counts +
               (xlim[1]+xlim[2]-(xlim[2]-xlim[1])*length(grp_info[[i]])/counts)/2
      result[grp_info[[i]],1] <- xgrid[2:(length(xgrid)-1)]
    } else {
      xgrid <- seq(0,xlim[2]-xlim[1],length.out=length(grp_info[[i]])+2)*tol +
               (xlim[1]+xlim[2]-(xlim[2]-xlim[1])*tol)/2
      result[grp_info[[i]],1] <- xgrid[2:(length(xgrid)-1)]
    }
  }
  grp_info[[n]] <- which(lambda>=breaks[n])
  if (length(grp_info[[n]])/counts>tol) {
    xgrid <- seq(0,xlim[2]-xlim[1],length.out=length(grp_info[[n]])+2)*length(grp_info[[n]])/counts +
             (xlim[1]+xlim[2]-(xlim[2]-xlim[1])*length(grp_info[[n]])/counts)/2
    result[grp_info[[n]],1] <- xgrid[2:(length(xgrid)-1)]
  } else {
    xgrid <- seq(0,xlim[2]-xlim[1],length.out=length(grp_info[[n]])+2)*tol +
      (xlim[1]+xlim[2]-(xlim[2]-xlim[1])*tol)/2
    result[grp_info[[n]],1] <- xgrid[2:(length(xgrid)-1)]
  }
  
  # optimizing horizontal axis by randomly swithing pairs
  if (ntrial>0) {
    for (j in seq(ntrial)) {
      random_group <- sample(seq(n),1)
      random_pair <- sample(grp_info[[random_group]],2)
      s1 <- result[random_pair[1],]
      s2 <- result[random_pair[2],]
      d1 <- c(s2[1],s1[2])
      d2 <- c(s1[1],s2[2])
      neighbor_01 <- as.vector(neighbors(graph,random_pair[1],mode="all"))
      neighbor_02 <- as.vector(neighbors(graph,random_pair[2],mode="all"))
      old_dist <- 0
      new_dist <- 0
      for (k1 in seq(length(neighbor_01))) {
        old_dist <- old_dist + sqrt(sum((s1-result[neighbor_01[k1],])^2))
        new_dist <- new_dist + sqrt(sum((d1-result[neighbor_01[k1],])^2))
      }
      for (k2 in seq(length(neighbor_02))) {
        old_dist <- old_dist + sqrt(sum((s2-result[neighbor_02[k2],])^2))
        new_dist <- new_dist + sqrt(sum((d2-result[neighbor_02[k2],])^2))
      }
      if (new_dist < old_dist) {
        result[random_pair[1],] <- d1
        result[random_pair[2],] <- d2
      }
    }
  }
  
  return(result)
}

# Function to generate customized graph layout based on fitted lambdas. Genes are divided into different 
# chunks based on percentiles of their lambdas. The vertical axis correspond to log lambdas, and the 
# horizontal axis within each chunk are evenly distributed and optimized (by randomly switch pairs of 
# horizontal coordinates within each chunks) so that the total plot length of edges are minimized.
#
# Input parameters:
#   graph: igraph object to be plotted
#   lambda: vector of fitted lambdas
#   xlim: 2-element vector that specify plot ranges for horizontal axis
#   ylim: 2-element vector that specify plot ranges for vertical axis
#   breaks: vector of percentiles to divide genes into different chunks
#   ntrials: number of trials to perform randomly switching of paris of horizontal coordinates
#   counts: a normalization parameter, together with `tol`, to scale the xlim within each chunk
#   tol: a normalization parameter, together with `counts`, to scale the xlim within each chunk
#
# Output:
#   result: 2-column matrix of coordinates as the plot layout
get_layout_loglambda <- function(graph,lambda,xlim,ylim,breaks,ntrial,counts,tol) {
  n <- length(breaks)
  result <- matrix(0,nrow=gorder(graph),ncol=2)
  loglambda <- sapply(lambda,function(x) log(max(min(x,1.5),0.03)))
  result[,2] <- loglambda/max(loglambda)*(ylim[2]-ylim[1])+ylim[1]
  grp_info <- vector("list",n)
  
  # assign initial horizontal coordinates within each chunk
  for (i in seq(n-1)) {
    grp_info[[i]] <- which(lambda>=breaks[i] & lambda<breaks[i+1])
    if (length(grp_info[[i]])/counts>tol) {
      xgrid <- seq(0,xlim[2]-xlim[1],length.out=length(grp_info[[i]])+2)*length(grp_info[[i]])/counts +
        (xlim[1]+xlim[2]-(xlim[2]-xlim[1])*length(grp_info[[i]])/counts)/2
      result[grp_info[[i]],1] <- xgrid[2:(length(xgrid)-1)]
    } else {
      xgrid <- seq(0,xlim[2]-xlim[1],length.out=length(grp_info[[i]])+2)*tol +
        (xlim[1]+xlim[2]-(xlim[2]-xlim[1])*tol)/2
      result[grp_info[[i]],1] <- xgrid[2:(length(xgrid)-1)]
    }
  }
  grp_info[[n]] <- which(lambda>=breaks[n])
  if (length(grp_info[[n]])/counts>tol) {
    xgrid <- seq(0,xlim[2]-xlim[1],length.out=length(grp_info[[n]])+2)*length(grp_info[[n]])/counts +
      (xlim[1]+xlim[2]-(xlim[2]-xlim[1])*length(grp_info[[n]])/counts)/2
    result[grp_info[[n]],1] <- xgrid[2:(length(xgrid)-1)]
  } else {
    xgrid <- seq(0,xlim[2]-xlim[1],length.out=length(grp_info[[n]])+2)*tol +
      (xlim[1]+xlim[2]-(xlim[2]-xlim[1])*tol)/2
    result[grp_info[[n]],1] <- xgrid[2:(length(xgrid)-1)]
  }
  
  # optimizing horizontal axis by randomly swithing pairs
  if (ntrial>0) {
    for (j in seq(ntrial)) {
      random_group <- sample(seq(n),1)
      random_pair <- sample(grp_info[[random_group]],2)
      s1 <- result[random_pair[1],]
      s2 <- result[random_pair[2],]
      d1 <- c(s2[1],s1[2])
      d2 <- c(s1[1],s2[2])
      neighbor_01 <- as.vector(neighbors(graph,random_pair[1],mode="all"))
      neighbor_02 <- as.vector(neighbors(graph,random_pair[2],mode="all"))
      old_dist <- 0
      new_dist <- 0
      for (k1 in seq(length(neighbor_01))) {
        old_dist <- old_dist + sqrt(sum((s1-result[neighbor_01[k1],])^2))
        new_dist <- new_dist + sqrt(sum((d1-result[neighbor_01[k1],])^2))
      }
      for (k2 in seq(length(neighbor_02))) {
        old_dist <- old_dist + sqrt(sum((s2-result[neighbor_02[k2],])^2))
        new_dist <- new_dist + sqrt(sum((d2-result[neighbor_02[k2],])^2))
      }
      if (new_dist < old_dist) {
        result[random_pair[1],] <- d1
        result[random_pair[2],] <- d2
      }
    }
  }
  
  return(result)
}


