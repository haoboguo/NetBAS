---
title: "HallMark Gene Sets Summary"
author: "HBG"
date: "12/7/2018"
output: pdf_document
---
  
#This script perform GO enrichment using NetBAS
#for 51 PDAC cell high BF genes
  
```{r}

rm(list=ls())
library(igraph)
library(gplots)

#list of 50 hallmark sets
list.file <- read.csv("list", header=F, stringsAsFactors=F)
hmk.set.list <- list.file$V1

#human PIN
human.pin <- read.csv("../../human.pin.csv", header=T, stringsAsFactors = F)
geneA <- human.pin$geneA
geneB <- human.pin$geneB

#All genes appear in PIN
all.list <- unique(c(geneA, geneB))

for (i in 1:length(hmk.set.list)) {
  ##the gene list and sub-network
  file.name <- paste("hmk.set.list[i]", "/", "hmk.set.list[i]", ".csv", sep="")
  geneset <- read.csv(file.name, header=T, stringsAsFactors=F)
  gele.list <- geneset$gene
  gene.number <- length(gene.list)   # this is the original number of genes the set
  #Constructing the sub network fromn the gene set
  subA <- geneA[which((geneA %in% gene.list) & (geneB %in% gene.list))]
  subB <- geneB[which((geneA %in% gene.list) & (geneB %in% gene.list))]
  genes.pin <- unique(c(subA, subB))
  gene.number.pin <- length(unique(genes.pin)) # this is the number of vertices of the sub-network
  sub.web <- data.frame(cbind(subA, subB))
  sub.edge.number <- length(sub.web[,1])  #this is the total edges of the sub-network
  
  # using igraph to build a graph of the sub-network
  sub.graph <- graph.data.frame(sub.web, directed = F)
  
  #number and degree of the largest clique
  largest.clique.number <- length(largest_cliques(sub.graph))
  largest.clique.degree <- length(largest_cliques(sub.graph)[[1]])
  
  #normalized values for betweenness, degree of vertices
  sub.betweenness <- betweenness(sub.graph, v=V(sub.graph), directed=F, weights=NULL, normalized=T)
  sub.degree <- degree(sub.graph, v=V(sub.graph), normalized = T)
  sub.closeness <- closeness(sub.graph, v=V(sub.graph), weights=NULL, normalized=T)
  number.of.clique <- clique.number(sub.graph) #total number of cliques
  sub.largest.clique <- largest_cliques(sub.graph)
  number.of.largest.clique <- length(sub.largest.clique)
  degree.of.largest.clique <- length(sub.largest.clique[[1]])
  cor.deg.bet <- cor(sub.degree, sub.betweenness) #correlation of degree/betweenness
  cor.deg.clo <- cor(sub.degree, sub.closeness) #correlation of degree/closeness
  cor.bet.clo <- cor(sub.betweenness, sub.closeness) #correlation of betweenness/closeness
  
  ##the enrichments result from NetBAS and DAVID
  
}

```