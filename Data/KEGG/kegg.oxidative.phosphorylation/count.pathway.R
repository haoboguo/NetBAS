rm(list=ls())
library(plyr)

network <- read.csv("/home/hguo/GitHub/NetBAS/Data/human.pin.csv", header=T, stringsAsFactors=F)
geneA <- network$geneA
geneB <- network$geneB

#all GO terms
pw.file <- read.csv("/home/hguo/GitHub/NetBAS/Data/human.pathway.csv", sep="\t", header=T, stringsAsFactors=F)
pw.gene <- pw.file$Gene
pw.pathway <- pw.file$Pathway
pw.source <- pw.file$Source

kegg.pw.cat <- unique(sort(pw.pathway[which(pw.source == "KEGG")]))
kegg.dim <- length(kegg.pw.cat)

gene.file <- read.csv("op.csv", header=T)
gene.list <- gene.file$gene 

##The background is regarded as the whole PIN
##for the rest of the list, the z-scores will be opposite
##to those for the gene set

glA <- geneB[which(geneA %in% gene.list)]
glB <- geneA[which(geneB %in% gene.list)]
gl.all <- c(glA, glB)
gl.count <- count(gl.all)
gl.orf <- gl.count$x
gl.freq <- gl.count$freq

# the gene list
kegg.vec <- numeric(length=kegg.dim)

    for (i in 1:length(gl.orf)) {
      orf.ith <- gl.orf[i]
      orf.freq <- gl.freq[i]
      orf.kegg.term <- pw.pathway[which((pw.gene %in% orf.ith) & (pw.source == "KEGG"))]
      for (j in 1:length(orf.kegg.term)) {
        na <- which(kegg.pw.cat %in% orf.kegg.term[j])
        kegg.vec[na] <- kegg.vec[na] + orf.freq
      }
  }

kegg.A <- matrix(kegg.vec, nrow = kegg.dim, ncol = 1)

write.table(kegg.A, file="HMK.kegg.csv", sep=",", col.names=F, row.names=F, quote=F)

