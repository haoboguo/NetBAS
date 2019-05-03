rm(list=ls())

human.pin <- read.csv("../../human.pin.csv", header=T, stringsAsFactors = F)
geneA <- human.pin$geneA
geneB <- human.pin$geneB

#All genes appear in PIN
all.list <- unique(c(geneA, geneB))

#the Gene list
geneset <- read.csv("HALLMARK_ESTROGEN_RESPONSE_LATE.csv", header=T, stringsAsFactors = F)
gene.list <- geneset$gene

#number of genes in the set
length(gene.list)

#number of genes of the set appear in the PIN
length(which(gene.list %in% all.list))

#Constructing the sub network fromn the gene set
subA <- geneA[which((geneA %in% gene.list) & (geneB %in% gene.list))]
subB <- geneB[which((geneA %in% gene.list) & (geneB %in% gene.list))]
sub.web <- data.frame(cbind(subA, subB))
length(sub.web[,1])

library(igraph)
library(gplots)
library(ggraph)
#'%ni%' <- Negate('%in%')
#connected.set <- unique(c(subA, subB))
#isolated <- gene.list[which(gene.list %ni% connected.set)]

sub.graph <- graph.data.frame(sub.web, directed = F)

dendro <- cluster_fast_greedy(sub.graph)
pdf("HALLMARK_ESTROGEN_RESPONSE_LATE.dendro.pdf", width=10, height=4, paper='special')
plot_dendrogram(dendro, rect=0, ann=F)
dev.off()

cluster <- cluster_walktrap(sub.graph)
pdf("HALLMARK_ESTROGEN_RESPONSE_LATE.custer.pdf", width=5, height=5, paper='special')
plot(cluster, sub.graph, modularity=T, vertex.label=NA, vertex.size=3)
dev.off()

