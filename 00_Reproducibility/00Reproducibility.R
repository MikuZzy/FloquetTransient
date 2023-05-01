# This file checks the reproducibility of the current data set by conducting differential cycling detection
# with another previously collected data set. Genes detected as significantly differential cycling under either
# temperatures are excluded from downstream analysis. 

setwd("~/R/novogene_trimmed")
library("limorhyde")
library("limma")
library("data.table")

#######################
# gene names
tpm2021 <- read.table(file = 'Data_novogene2021_trimmed_tpm.tsv', header = TRUE)
gene_keys <- read.csv(file = 'Data_SampleKey_v3.csv', header = TRUE)
tpmFB <- read.table(file = 'Data_fatbody2012_trimmed_rsem_tpm.tsv', header = TRUE)
setwd("~/R/novogene_trimmed/00_Reproducibility")

gene_suffix <- data.frame(
  "geneid" = sort(intersect(tpm2021$geneid,tpmFB$gene))
)
gene_suffix$genename <- sapply(gene_suffix$geneid,
                               function(x,y) y$gene_name[which(y$geneid==x)],
                               y=tpm2021)
gene_suffix$gene <- apply(gene_suffix,1,function(x) paste(x[1],x[2],sep="_"))
coreclock <- c("FBgn0023076_Clk","FBgn0016076_vri","FBgn0003068_per","FBgn0014396_tim")

#######################
# V2 data, current data set
v2_18C_day5 <- which(gene_keys$Temperature=="18C" & gene_keys$Day==5)
v2_25C_day5 <- which(gene_keys$Temperature=="25C" & gene_keys$Day==5)

metadata_v2_18C <- data.frame(
  "sample" = gene_keys$Label[v2_18C_day5],
  "cond" = "V2_18C",
  "time" = gene_keys$Hour[v2_18C_day5]
)
metadata_v2_25C <- data.frame(
  "sample" = gene_keys$Label[v2_25C_day5],
  "cond" = "V2_25C",
  "time" = gene_keys$Hour[v2_25C_day5]
)

ind_v2 <- as.vector(unlist(sapply(gene_suffix$geneid,
                                  function(x,y) which(y==x), y=tpm2021$geneid)))
data_v2_18C <- as.matrix(tpm2021[ind_v2,v2_18C_day5+2])
rownames(data_v2_18C) <- gene_suffix$gene
data_v2_25C <- as.matrix(tpm2021[ind_v2,v2_25C_day5+2])
rownames(data_v2_25C) <- gene_suffix$gene

#######################
# V1 data, a previously collected data set
v1_18C_day5 <- seq(2,25,1)
v1_25C_day5 <- seq(26,49,1)

metadata_v1_18C <- data.frame(
  "sample" = colnames(tpmFB)[v1_18C_day5],
  "cond" = "V1_18C",
  "time" = rep(seq(2,24,2),2)
)
metadata_v1_25C <- data.frame(
  "sample" = colnames(tpmFB)[v1_25C_day5],
  "cond" = "V1_25C",
  "time" = rep(seq(2,24,2),2)
)

ind_v1 <- as.vector(unlist(sapply(gene_suffix$geneid,
                                  function(x,y) which(y==x), y=tpmFB$gene)))
data_v1_18C <- as.matrix(tpmFB[ind_v1,v1_18C_day5])
rownames(data_v1_18C) <- gene_suffix$gene
data_v1_25C <- as.matrix(tpmFB[ind_v1,v1_25C_day5])
rownames(data_v1_25C) <- gene_suffix$gene

#######################
# genes with median >5 in either temperature in v2 dataset, and z-score transforming
ind_filter_18C <- as.vector(apply(data_v2_18C,1,function(x) median(x)>5))
ind_filter_25C <- as.vector(apply(data_v2_25C,1,function(x) median(x)>5))
ind_filter <- (ind_filter_18C+ind_filter_25C)>0

data_v2_18C <- t(apply(data_v2_18C[ind_filter,],1,function(x) (x-mean(x))/sd(x)*sqrt((length(x)-1)/length(x))))
data_v2_25C <- t(apply(data_v2_25C[ind_filter,],1,function(x) (x-mean(x))/sd(x)*sqrt((length(x)-1)/length(x))))
data_v1_18C <- t(apply(data_v1_18C[ind_filter,],1,function(x) (x-mean(x))/sd(x)*sqrt((length(x)-1)/length(x))))
data_v1_25C <- t(apply(data_v1_25C[ind_filter,],1,function(x) (x-mean(x))/sd(x)*sqrt((length(x)-1)/length(x))))

#######################
# limorhyde dataset and metadata
data_18C <- cbind(data_v1_18C,data_v2_18C)
metadata_18C <- rbind.data.frame(metadata_v1_18C,metadata_v2_18C)
data_25C <- cbind(data_v1_25C,data_v2_25C)
metadata_25C <- rbind.data.frame(metadata_v1_25C,metadata_v2_25C)

## 18C fit
metadata_18C <- cbind(metadata_18C, limorhyde(metadata_18C$time, 'time_'))
design_18C <- model.matrix(~ cond*(time_cos+time_sin), data=metadata_18C)
fit_18C <- lmFit(data_18C, design_18C)
fit_18C <- eBayes(fit_18C, trend=TRUE)
drLimma_18C <- data.table(topTable(fit_18C,coef=5:6,number=Inf),keep.rownames=TRUE)
setnames(drLimma_18C,'rn','gene_id')
setorderv(drLimma_18C,'P.Value')
drLimma_18C[which(drLimma_18C$gene_id %in% coreclock),]
drLimma_18C[which(drLimma_18C$P.Value<=0.05),]

## 25C fit
metadata_25C <- cbind(metadata_25C, limorhyde(metadata_25C$time, 'time_'))
design_25C <- model.matrix(~ cond*(time_cos+time_sin), data=metadata_25C)
fit_25C <- lmFit(data_25C, design_25C)
fit_25C <- eBayes(fit_25C, trend=TRUE)
drLimma_25C <- data.table(topTable(fit_25C,coef=5:6,number=Inf),keep.rownames=TRUE)
setnames(drLimma_25C,'rn','gene_id')
setorderv(drLimma_25C,'P.Value')
drLimma_25C[which(drLimma_25C$gene_id %in% coreclock),]

#######################
# not significantly differentially cycled genes
alpha <- 0.05
geneid_ndc <- sort(intersect(drLimma_18C$gene_id[which(drLimma_18C$P.Value>alpha)],
                             drLimma_25C$gene_id[which(drLimma_25C$P.Value>alpha)]))
result <- data.frame(
  "geneid" = as.vector(sapply(geneid_ndc,function(x) strsplit(x,split="_")[[1]][1])),
  "gene_name" = as.vector(sapply(geneid_ndc,function(x) strsplit(x,split="_")[[1]][2]))
)
write.csv(result,file="00ReproducibleGenes.csv")

