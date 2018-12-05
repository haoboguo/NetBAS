rm(list=ls())

human.pin <- read.csv("../../human.pin.csv", header=T, stringsAsFactors = F)
geneA <- human.pin$geneA
geneB <- human.pin$geneB

#All genes appear in PIN
all.list <- unique(c(geneA, geneB))

#the Gene list
panc.file <- read.csv("../../../rnf43.csv",header=TRUE,stringsAsFactors=F)
panc.gene <- panc.file$gene
panc.z <- panc.file$zPanc_vs_other
gene.list <- panc.gene[which(panc.z > 2)]

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

sub.graph <- graph.data.frame(sub.web, directed = F)

pdf("PANC17.full.pdf")
plot(sub.graph, vertex.label=NA, vertex.size=6)
dev.off()

#the largest clique of the sub network
sub.largest.clique <- largest_cliques(sub.graph)

#total number of largest clques
length(sub.largest.clique)

#number of genes in the largest cliques
sub.clique.1 <- sub.largest.clique[[1]]
degree <- length(sub.clique.1)
degree

#clique 1
sub.c1.graph <- graph.full(degree)
V(sub.c1.graph)$name <- V(sub.graph)$name[sub.clique.1]

pdf("PANC17.clique.1.pdf")
plot(sub.c1.graph)
dev.off()
