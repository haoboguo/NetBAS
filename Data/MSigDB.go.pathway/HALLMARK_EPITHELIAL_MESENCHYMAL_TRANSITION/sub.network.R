rm(list=ls())

human.pin <- read.csv("../../human.pin.csv", header=T, stringsAsFactors = F)
geneA <- human.pin$geneA
geneB <- human.pin$geneB

#All genes appear in PIN
all.list <- unique(c(geneA, geneB))

#the Gene list
geneset <- read.csv("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.csv", header=T, stringsAsFactors = F)
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
'%ni%' <- Negate('%in%')
connected.set <- unique(c(subA, subB))
isolated <- gene.list[which(gene.list %ni% connected.set)]

sub.graph <- graph.data.frame(sub.web, directed = F)
sub.graph.iso <- sub.graph %>% add_vertices(length(isolated))

pdf("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.full.pdf")
plot(sub.graph.iso, vertex.label=NA, vertex.size=6)
dev.off()

#the largest clique of the sub network
sub.largest.clique <- largest_cliques(sub.graph.iso)

#total number of largest clques
length(sub.largest.clique)

#number of genes in the largest cliques
sub.clique.1 <- sub.largest.clique[[1]]
degree <- length(sub.clique.1)
degree

#clique 1
sub.c1.graph <- graph.full(degree)
V(sub.c1.graph)$name <- V(sub.graph)$name[sub.clique.1]

pdf("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.clique1.pdf")
plot(sub.c1.graph)
dev.off()
