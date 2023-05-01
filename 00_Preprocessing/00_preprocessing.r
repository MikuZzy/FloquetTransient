# This file preprocesses the tpm data with the same output format.

setwd("~/R/novogene_trimmed/00_Preprocessing")

novogene_Ntrim_dedup <- read.table(file='Data_novogene2021_dedup_rsem_tpm.tsv',header=TRUE)
novogene_trim_dedup <- read.table(file='novogene2021_trimmed_dedup_tpm.tsv',header=TRUE)
novogene_trim_Ndedup <- read.table(file='novogene2021_trimmed_tpm.tsv',header=TRUE)
novogene_isoform <- read.table(file='novogene2021_trimmed_dedup_isoforms_tpm.tsv',header=TRUE)

DR_Ntrim <- read.table(file='Data_fatbody_genes_rsem_tpm.tsv',header=TRUE)
DR_trim <- read.table(file='fatbody2012_trimmed_rsem_tpm.tsv',header=TRUE)

# Check that data files contain same collection of genes.

prod(novogene_trim_dedup$gene==novogene_trim_Ndedup$gene)
prod(novogene_trim_dedup$gene==DR_trim$gene)

gene_suffix01 <- data.frame(
  'geneid' = as.vector(sapply(novogene_trim_dedup$gene,
                              function(x) strsplit(x,split='__')[[1]][1])),
  'gene_name' = as.vector(sapply(novogene_trim_dedup$gene,
                                 function(x) strsplit(x,split='__')[[1]][2]))
)
gene_suffix02 <- data.frame(
  'gene' = as.vector(sapply(DR_trim$gene,
                            function(x) strsplit(x,split='__')[[1]][1]))
)
gene_suffix03 <- data.frame(
  'transcriptid' = as.vector(sapply(novogene_isoform$transcript,
                                    function(x) strsplit(x,split='__')[[1]][1])),
  'geneid' = as.vector(sapply(novogene_isoform$transcript,
                              function(x) strsplit(x,split='__')[[1]][2])),
  'gene_name' = as.vector(sapply(novogene_isoform$transcript,
                                 function(x) strsplit(x,split='__')[[1]][3]))
)

novogene_trim_dedup <- cbind.data.frame(gene_suffix01,
                                        novogene_trim_dedup[,2:ncol(novogene_trim_dedup)])
novogene_trim_Ndedup <- cbind.data.frame(gene_suffix01,
                                         novogene_trim_Ndedup[,2:ncol(novogene_trim_Ndedup)])
DR_trim <- cbind.data.frame(gene_suffix02,DR_trim[,2:ncol(DR_trim)])
novogene_isoform <- cbind.data.frame(gene_suffix03,
                                     novogene_isoform[,2:ncol(novogene_isoform)])

# Check data frame format after preprocessing.

prod(colnames(novogene_trim_dedup)==colnames(novogene_Ntrim_dedup))
prod(colnames(novogene_trim_Ndedup)==colnames(novogene_Ntrim_dedup))
prod(colnames(DR_trim)==colnames(DR_Ntrim))

setwd("~/R/novogene_trimmed")
write.table(novogene_trim_dedup,file='Data_novogene2021_trimmed_dedup_tpm.tsv',
            row.names=FALSE,sep="\t",quote=FALSE)
write.table(novogene_trim_Ndedup,file='Data_novogene2021_trimmed_tpm.tsv',
            row.names=FALSE,sep="\t",quote=FALSE)
write.table(DR_trim,file='Data_fatbody2012_trimmed_rsem_tpm.tsv',
            row.names=FALSE,sep="\t",quote=FALSE)
write.table(novogene_isoform,file='Data_novogene2021_trimmed_dedup_isoforms_tpm.tsv',
            row.names=FALSE,sep="\t",quote=FALSE)

