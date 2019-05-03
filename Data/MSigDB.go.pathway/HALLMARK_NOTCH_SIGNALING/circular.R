rm(list=ls())

human.pin <- read.csv("../../human.pin.csv", header=T, stringsAsFactors = F)
geneA <- human.pin$geneA
geneB <- human.pin$geneB

#All genes appear in PIN
all.list <- unique(c(geneA, geneB))

#the Gene list
geneset <- read.csv("HALLMARK_NOTCH_SIGNALING.csv", header=T, stringsAsFactors = F)
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
'%ni%' <- Negate('%in%')
connected.set <- unique(c(subA, subB))
isolated <- gene.list[which(gene.list %ni% connected.set)]

sub.graph <- graph.data.frame(sub.web, directed = F)
sub.graph.iso <- sub.graph %>% add_vertices(length(isolated))

pdf("HALLMARK_NOTCH_SIGNALING.circular.pdf", width=8, height=7.5, paper='special')
ggraph(sub.graph.iso, layout='linear', circular=T) + geom_edge_arc() + geom_node_point() +
      theme_graph(foreground = 'steelblue', fg_text_colour = 'white')
dev.off()

