# This file performs network analysis of concordant genes. It outputs figures of various plots as well as a table
# summary of each pathway.

setwd("~/R/novogene_trimmed")
results <- read.csv(file='01Result_concordant_details_repro.csv',header=TRUE)
results_fullpara <- read.csv(file='01Result_filtered_fullparafit.csv',header=TRUE)

setwd("~/R/novogene_trimmed/05_GRN")
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.0",.libPaths()))
library('graphite')
library('igraph')
library('colorspace')
source("nwk_srcfuncs.r")

# Get pathway from reactome.
d <- pathways("dmelanogaster","reactome")
num_pathway <- length(names(d))

GRN_overview <- matrix(NA,nrow=num_pathway,ncol=9)
colnames(GRN_overview) <- c("Name","NumEdge_nobinding","NumVert_nobinding",
                            "NumEdge_w/lambda","NumVert_w/lambda",
                            "NumEdge_Induce","NumVert_Induce",
                            "VertFrac_w/lambda","VertFrac_converage")

GRN_overview <- as.data.frame(GRN_overview)
GRN_overview$Name <- names(d)
GRN_graph_data <- data.frame(matrix(ncol=9,nrow=0))
colnames(GRN_graph_data) <- c("src_type","src","dest_type","dest","direction","type","src_lambda","dest_lambda","induce")
tol <- 0.1
lambda_dist <- rep(NA,240000)
index_towrite <- 1

# Loop through all pathways. Select concordant vertices, and edges that are not 'binding' and whose lambda 
# differences are not smaller than tol = -0.1.
for (i in seq(num_pathway)) {
  p <- d[[names(d)[i]]]
  p <- convertIdentifiers(p,"FLYBASE")
  e_info <- graphite::edges(p)
  ind_notbinding <- which(e_info$type!="Binding")
  e_info <- e_info[ind_notbinding,]
  e_info <- unique(e_info)
  
  # concordant vertices
  start <- as.numeric(sapply(e_info$src,is.element,results$geneid))
  end <- as.numeric(sapply(e_info$dest,is.element,results$geneid))
  e_filtered <- e_info[which(start*end==1),]
  
  # edges with lambda differences no smaller than tol = -0.1
  e_filtered$src_lambda <- as.vector(sapply(e_filtered$src,lambda_from_geneid,results=results))
  e_filtered$dest_lambda <- as.vector(sapply(e_filtered$dest,lambda_from_geneid,results=results))
  if (nrow(e_filtered) > 0) {
    e_filtered$induce <- as.numeric(e_filtered$src_lambda>=e_filtered$dest_lambda-tol)
  }
  e_filtered$src <- as.vector(sapply(e_filtered$src,genename_from_geneid,results=results))
  e_filtered$dest <- as.vector(sapply(e_filtered$dest,genename_from_geneid,results=results))
  
  e_induce <- e_filtered[which(e_filtered$induce==1),]
  
  # write output
  GRN_overview$NumEdge_nobinding[i] <- nrow(e_info)
  GRN_overview$NumVert_nobinding[i] <- length(unique(c(e_info$src,e_info$dest)))
  
  GRN_overview$`NumEdge_w/lambda`[i] <- nrow(e_filtered)
  GRN_overview$`NumVert_w/lambda`[i] <- length(unique(c(e_filtered$src,e_filtered$dest)))
  if (GRN_overview$NumVert_nobinding[i]!=0) {
    GRN_overview$`VertFrac_w/lambda`[i] <- GRN_overview$`NumVert_w/lambda`[i] / GRN_overview$NumVert_nobinding[i]
  } else {
    GRN_overview$`VertFrac_w/lambda`[i] <- 0
  }
  
  GRN_overview$NumEdge_Induce[i] <- nrow(e_induce)
  GRN_overview$NumVert_Induce[i] <- length(unique(c(e_induce$src,e_induce$dest)))
  if (GRN_overview$`NumVert_w/lambda`[i]!=0) {
    GRN_overview$VertFrac_converage[i] <- GRN_overview$NumVert_Induce[i] / GRN_overview$`NumVert_w/lambda`[i]
  } else {
    GRN_overview$VertFrac_converage[i] <- 0
  }
  if (nrow(e_filtered) > 0) {
    lambda_dist[index_towrite:(index_towrite+nrow(e_filtered)-1)] <- e_filtered$src_lambda - e_filtered$dest_lambda
    index_towrite <- index_towrite + nrow(e_filtered)
  }
  if (nrow(e_induce)>0) {
    GRN_graph_data <- rbind.data.frame(GRN_graph_data,e_induce)
  }
  message(i,"/",num_pathway,', numVert = ',length(unique(c(GRN_graph_data$src,GRN_graph_data$dest))),
          ', numEdge = ',nrow(GRN_graph_data),', write to index ',index_towrite)
}

# table summary of pathways
GRN_overview <- GRN_overview[order(GRN_overview$VertFrac_converage,
                                   GRN_overview$`VertFrac_w/lambda`,
                                   GRN_overview$NumEdge_nobinding,
                                   decreasing=TRUE),]
write.csv(GRN_overview,file="05_GRN_overview_tol.csv")
save.image(file="05GRN_tol.Rdata")

# plot distributions of lambda differences of all non-binding edges 
lambda_dist[is.nan(lambda_dist)] <- 0
lambda_dist <- as.vector(na.omit(lambda_dist))
lambda_dist_plot <- sapply(lambda_dist, function(x) max(min(x,1.5),-1.5))
mean(lambda_dist>=-tol)

h1 <- hist(lambda_dist_plot,breaks=seq(-1.5,1.5,0.1), xlab='difference of lambda through directed edges', 
           main="bar width = 0.10")
h2 <- hist(lambda_dist_plot,breaks=seq(-1.5,1.5,0.05), xlab='difference of lambda through directed edges', 
           main="bar width = 0.05")
h3 <- hist(lambda_dist_plot,breaks=seq(-1.5,1.5,0.01), xlab='difference of lambda through directed edges', 
           main="bar width = 0.01")

GRN_graph <- get_graph(GRN_graph_data)
GRN_graph <- simplify(GRN_graph,remove.loops=FALSE)
out_deg <- degree(GRN_graph,v=V(GRN_graph),mode="out",loops=TRUE,normalized=FALSE)
Vcolor <- as.vector(ifelse(out_deg>=10,"red",ifelse(out_deg>=5,"orange","yellow")))

######### Color by outer degrees, whole graph
dev.new(width=8, height=15, unit="in")
plot.new()
colours = c("red","orange","yellow")
labels = paste("label", 1:length(colours))
legend("center",legend=c(">=10","5~9","<=4"), col=colours, pch=19,cex=1.7,
       title="out degrees")

set.seed(1)
tkplot(GRN_graph,
       canvas.width = 1600,
       canvas.height = 900,
       vertex.color = Vcolor,
       vertex.size = 20,
       vertex.label.font = 2,
       vertex.label.cex = 1.5,
       edge.width	= 4,
       edge.arrow.size = 2,
       edge.color = Ecolor,
       layout = layout_nicely
)

write.csv(as.vector(sapply(V(GRN_graph)$name,geneid_from_genename,results=results)),
          file='05_whole_tol.csv')
save.image(file="05GRN_tol.Rdata")

######### subgraph with connections to core clock genes
vnames <- c()
for (i in V(GRN_graph)$name) {
  a1 <- distances(GRN_graph,v=i,to="per",mode=c("out"))
  a2 <- distances(GRN_graph,v=i,to="tim",mode=c("out"))
  a3 <- distances(GRN_graph,v=i,to="Clk",mode=c("out"))
  a4 <- distances(GRN_graph,v="per",to=i,mode=c("out"))
  a5 <- distances(GRN_graph,v="tim",to=i,mode=c("out"))
  a6 <- distances(GRN_graph,v="Clk",to=i,mode=c("out"))
  if (a1<Inf|a2<Inf|a3<Inf|a4<Inf|a5<Inf|a6<Inf) {
    vnames <- append(vnames,i)
  }
  message(which(V(GRN_graph)$name==i),'/',length(V(GRN_graph)),', nvert = ',length(vnames))
}

GRN_coreclock <- induced_subgraph(GRN_graph,vids=vnames)
color_coreclock_by_lambda <- as.vector(sapply(V(GRN_coreclock)$name,lambda_from_genename,results=results))
lambda_palette <- c("#FFFF00","#FFE600","#FFCC00","#FFB300","#FF9900","#FF8000","#FF6600",
                    "#FF4D00","#FF3300","#FF1A00","#FF0000")
Vcolor_coreclock_lambda <- unlist(lapply(color_coreclock_by_lambda,function(x,y) y[min(max(ceiling(x/0.1),1),11)],y=lambda_palette))
cycling_18_coreclock <- as.vector(sapply(V(GRN_coreclock)$name,cycling_from_genename,results=results,temp=18))
cycling_25_coreclock <- as.vector(sapply(V(GRN_coreclock)$name,cycling_from_genename,results=results,temp=25))
Vcolor_temp18_coreclock <- as.vector(ifelse(cycling_18_coreclock,"green","grey"))
Vcolor_temp25_coreclock <- as.vector(ifelse(cycling_25_coreclock,"green","grey"))
Ecolor_coreclock <- as.vector(ifelse(which_mutual(GRN_coreclock),"gray","black"))

# get customized layouts, vertical coordinates either by rank or log of lambdas
# set.seed(1)
# coord_coreclock <- get_layout(GRN_coreclock,lambda=color_coreclock_by_lambda,xlim=c(0,1),ylim=c(0,1),
#                               breaks=as.numeric(quantile(color_coreclock_by_lambda,c(0,0.5))),
#                               ntrial=10000,counts=10,tol=0.2)
set.seed(1)
coord_coreclock_loglambda <- get_layout_loglambda(GRN_coreclock,lambda=color_coreclock_by_lambda,xlim=c(0,1),
                                                  ylim=c(0,1),breaks=as.numeric(quantile(color_coreclock_by_lambda,c(0,0.5))),
                                                  ntrial=10000,counts=20,tol=0.4)
label_direction <- c(pi/2,-pi/2,-pi/2,-pi/2,-pi/2,pi/2,-pi/2,-pi/2,pi/2,
                     -pi/2,-pi/2,pi/2,-pi/2,pi/2,-pi/2,pi/2,-pi/2,-pi/2,
                     -pi/2,pi/2,-pi/2,-pi/2,-pi/2,-pi/2,pi/2,pi/2,-pi/2,
                     pi/2,-pi/2,-pi/2,pi/2,pi/2)

# manual adjustment of horizontal coordinates for better visulization
coord_coreclock_loglambda_mod <- coord_coreclock_loglambda
coord_coreclock_loglambda_mod[2,1] <- 0.33 #Rpn6
coord_coreclock_loglambda_mod[3,1] <- 0.19 #Prosalpha4
coord_coreclock_loglambda_mod[9,1] <- 0.40 #Rpt3
coord_coreclock_loglambda_mod[10,1] <- 0.35 #Prosbeta7
coord_coreclock_loglambda_mod[5,1] <- 0.25 #Rpn8
coord_coreclock_loglambda_mod[20,1] <- 0.3 #Usp14

coord_coreclock_loglambda_mod[11,1] <- 0.7 #Prosbeta3
coord_coreclock_loglambda_mod[8,1] <- 0.59 #Rpt5
coord_coreclock_loglambda_mod[6,1] <- 0.64 #Rpn9
coord_coreclock_loglambda_mod[4,1] <- 0.54 #Rpn3

# plots of legends and subgraphs
dev.new(width=8, height=15, unit="in")
plot.new()
colours = rev(lambda_palette)
labels = paste("label", 1:length(colours))
legend("center",legend=c(">1.0","0.9~1.0","0.8~0.9","0.7~0.8","0.6~0.7","0.5~0.6",
                         "0.4~0.5","0.3~0.4","0.2~0.3","0.1~0.2","0~0.1"), col=colours, pch=19,cex=1.7,
       title="lambda")

dev.new(width=8, height=10, unit="in")
plot.new()
legend("center",legend=c("Uni-directional Edges","Bi-directional Edges"), col=c("black","grey"),lty=1,lwd=4,
       cex=1.7, title="Edges")

tkplot(GRN_coreclock,
       canvas.width = 1200,
       canvas.height = 675,
       vertex.color = Vcolor_coreclock_lambda,
       vertex.size = 20,
       vertex.label.font = 2,
       vertex.label.cex = 2.2,
       vertex.label.dist = 1.2,
       vertex.label.color = "dodgerblue",
       vertex.label.degree = label_direction,
       edge.width	= 3,
       edge.arrow.size = 2,
       edge.color = Ecolor_coreclock,
       layout = coord_coreclock_loglambda_mod,
       margin = rep(0.2,4)
) 

dev.new(width=8, height=15, unit="in")
plot.new()
colours = c("green","grey")
labels = paste("label", 1:length(colours))
legend("center",legend=c("significantly cycling",
                         "not significantly cycling"), col=colours, pch=19,cex=1.7,
       title="cycling under 18C")

tkplot(GRN_coreclock,
       canvas.width = 1200,
       canvas.height = 675,
       vertex.color = Vcolor_temp18_coreclock,
       vertex.size = 20,
       vertex.label.font = 2,
       vertex.label.cex = 2.2,
       vertex.label.dist = 1.2,
       vertex.label.color = "dodgerblue",
       vertex.label.degree = label_direction,
       edge.width	= 3,
       edge.arrow.size = 2,
       edge.color = Ecolor_coreclock,
       layout = coord_coreclock_loglambda_mod,
       margin = rep(0.2,4)
) 

dev.new(width=8, height=10, unit="in")
plot.new()
colours = c("green","grey")
labels = paste("label", 1:length(colours))
legend("center",legend=c("significantly cycling",
                         "not significantly cycling"), col=colours, pch=19,cex=1.7,
       title="cycling under 25C")

tkplot(GRN_coreclock,
       canvas.width = 1200,
       canvas.height = 675,
       vertex.color = Vcolor_temp25_coreclock,
       vertex.size = 20,
       vertex.label.font = 2,
       vertex.label.cex = 2.2,
       vertex.label.dist = 1.2,
       vertex.label.color = "dodgerblue",
       vertex.label.degree = label_direction,
       edge.width	= 3,
       edge.arrow.size = 2,
       edge.color = Ecolor_coreclock,
       layout = coord_coreclock_loglambda_mod,
       margin = rep(0.2,4)
) 

write.csv(as.vector(sapply(V(GRN_coreclock)$name,geneid_from_genename,results=results)),
          file='05_coreclock_tol.csv')
save.image(file="05GRN_tol.Rdata")

################################
# phase of upstream and downstream genes connected to Clk
dist_out <- distances(GRN_coreclock, v="Clk", to = V(GRN_coreclock), mode="out")
dist_in <- distances(GRN_coreclock, v="Clk", to = V(GRN_coreclock), mode="in")
ind_connected <- which(dist_out!=Inf | dist_in!=Inf)
dist_out <- dist_out[ind_connected]
dist_in <- dist_in[ind_connected]

phase_data <- data.frame(
  "genename" = V(GRN_coreclock)$name[ind_connected],
  "distance" = dist_out,
  "cycle18C" = cycling_18_coreclock[ind_connected],
  "cycle25C" = cycling_25_coreclock[ind_connected],
  "phase18C" = rep(NA,length(ind_connected)),
  "phase25C" = rep(NA,length(ind_connected))
)
phase_data$distance[which(dist_out==Inf)] <- -dist_in[which(dist_out==Inf)]
phase_data <- phase_data[order(phase_data$distance),]
for (i in seq(nrow(phase_data))) {
  j <- which(results_fullpara$gene_name==phase_data$genename[i])
  if (phase_data$cycle18C[i]) {
    phase_data$phase18C[i] <- atan2(results_fullpara$sin_18C[j]+results_fullpara$sin_25C[j],
                                    results_fullpara$cos_18C[j]+results_fullpara$cos_25C[j])/pi*12
  }
  if (phase_data$cycle25C[i]) {
    phase_data$phase25C[i] <- atan2(results_fullpara$sin_25C[j],
                                    results_fullpara$cos_25C[j])/pi*12
  }
  message(i)
}

phase_data$distance_mod <- phase_data$distance
phase_data$distance_mod[which(phase_data$distance_mod>2)] <- "3 or more"
phase_data$distance_mod[which(phase_data$distance_mod==-1)] <- "-1 (pointing to Clk)"

phase_data_cycle18C <- phase_data[phase_data$cycle18C,c(1,2,3,5,7)]
phase_data_cycle18C <- phase_data_cycle18C[order(phase_data_cycle18C$distance),]
phase_data_cycle18C_projphase <- phase_data_cycle18C
phase_data_cycle18C_projphase$phase18C <-phase_data_cycle18C_projphase$phase18C+24
phase_data_cycle18C <- rbind.data.frame(phase_data_cycle18C,phase_data_cycle18C_projphase)
phase_data_cycle18C$temp <- '18C'
colnames(phase_data_cycle18C) <- c('genename','distance','cycle','phase','distance_mod','temp')

phase_data_cycle25C <- phase_data[phase_data$cycle25C,c(1,2,4,6,7)]
phase_data_cycle25C <- phase_data_cycle25C[order(phase_data_cycle25C$distance),]
phase_data_cycle25C_projphase <- phase_data_cycle25C
phase_data_cycle25C_projphase$phase25C <-phase_data_cycle25C_projphase$phase25C+24
phase_data_cycle25C <- rbind.data.frame(phase_data_cycle25C,phase_data_cycle25C_projphase)
phase_data_cycle25C$temp <- '25C'
colnames(phase_data_cycle25C) <- c('genename','distance','cycle','phase','distance_mod','temp')

library(ggplot2)
phase_data_plot <- rbind.data.frame(phase_data_cycle18C, phase_data_cycle25C)
phase_data_plot <- phase_data_plot[order(phase_data_plot$distance,phase_data_plot$temp),]
plot_order <- V(GRN_coreclock)$name[order(color_coreclock_by_lambda,decreasing=FALSE)]
plot_order <- plot_order[which(plot_order%in%phase_data_plot$genename)]

# genes ordered by rates
g3 <- ggplot(data = phase_data_plot) + 
  geom_point(aes(y=phase,x=genename,col=factor(distance_mod),fill=factor(distance_mod),shape=temp),size=5) +
  coord_flip() +  scale_x_discrete(limits=plot_order) +
  scale_y_continuous(breaks = seq(-12,36,by=6)) + xlab("Gene Names") + ylab("Fitted Phase under 25C and 18C") +
  scale_shape_manual(values=c('25C'=2,'18C'=25),breaks=c("25C","18C")) +
  guides(col=guide_legend(title="Distance to Clk"),shape=guide_legend(title="Temperature"),
         fill=guide_legend(title="Distance to Clk")) +
  expand_limits(y=c(-12,36))

save.image(file="05GRN_tol.Rdata")

