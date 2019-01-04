rm(list=ls())
#library(plyr)
library(gplots)
library(microbenchmark)
library(matrixStats)

#KEGG pathways
pw.file <- read.csv("/home/hguo/GitHub/NetBAS/Data/human.pathway.csv", sep="\t", header=T, stringsAsFactors=F)
pw.gene <- pw.file$Gene
pw.pathway <- pw.file$Pathway
pw.source <- pw.file$Source

kegg.pw.cat <- unique(sort(pw.pathway[which(pw.source == "KEGG")]))
kegg.dim <- length(kegg.pw.cat)

kegg.hspin <- matrix(as.numeric(unlist(read.table("HMK.kegg.csv", header=F, sep=","))), ncol=1)
kegg.obs <- c(kegg.hspin)

kegg.perm <- c()
for (i in 1:1000) {
    kegg.name <- paste("ms02.", i, ".kegg.csv", sep="")
    kegg.mat <- matrix(as.numeric(unlist(read.table(kegg.name, header=F, sep=","))), ncol=1)
    kegg.perm <- rbind(kegg.perm, c(kegg.mat))
}

kegg.mean <- colMeans(kegg.perm)
kegg.std <- colSds(kegg.perm)

kegg.zscore <- round((kegg.obs - kegg.mean)/kegg.std, 3)

kegg.order <- order(-kegg.zscore)

z.kegg <- matrix(kegg.zscore[kegg.order], nrow=kegg.dim)

rownames(z.kegg) <- kegg.pw.cat[kegg.order]

kegg.enrich.list <- rownames(z.kegg)

kegg.enriched.terms <- c("Pathway", "Z-score")
for (i in 1:length(kegg.enrich.list)) {
  id <- as.character(kegg.enrich.list[i])
  z.gene <- z.kegg[i]
  kegg.enriched.terms <- rbind(kegg.enriched.terms, c(id, z.gene))
}
kegg.enriched.terms
write.table(kegg.enriched.terms, file="Kegg-OP.kegg.enriched.txt", sep="\t", row.names = F, col.names = F, quote=T)

