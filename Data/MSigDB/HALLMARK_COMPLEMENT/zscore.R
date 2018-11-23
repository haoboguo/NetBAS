rm(list=ls())
library(plyr)
library(gplots)
library(microbenchmark)
library(matrixStats)
library(GO.db)

## the GO terms for biological processes (BP)
bp.GOcategory.file <- read.csv("/home/hguo/GitHub/NetBAS/Data/human.bp.term.csv",header=TRUE, stringsAsFactors=F)
bp.go.cat <- bp.GOcategory.file$GO.term
bp.dim <- length(bp.go.cat)
cc.GOcategory.file <- read.csv("/home/hguo/GitHub/NetBAS/Data/human.cc.term.csv",header=TRUE, stringsAsFactors=F)
cc.go.cat <- cc.GOcategory.file$GO.term
cc.dim <- length(cc.go.cat)
mf.GOcategory.file <- read.csv("/home/hguo/GitHub/NetBAS/Data/human.mf.term.csv",header=TRUE, stringsAsFactors=F)
mf.go.cat <- mf.GOcategory.file$GO.term
mf.dim <- length(mf.go.cat)

bp.hspin <- matrix(as.numeric(unlist(read.table("HMK.bp.csv", header=F, sep=","))), nrow=bp.dim, ncol=2)
bp.obs <- c(bp.hspin)
cc.hspin <- matrix(as.numeric(unlist(read.table("HMK.cc.csv", header=F, sep=","))), nrow=cc.dim, ncol=2)
cc.obs <- c(cc.hspin)
mf.hspin <- matrix(as.numeric(unlist(read.table("HMK.mf.csv", header=F, sep=","))), nrow=mf.dim, ncol=2)
mf.obs <- c(mf.hspin)

bp.perm <- c()
cc.perm <- c()
mf.perm <- c()
for (i in 1:1000) {
    bp.name <- paste("ms02.", i, ".bp.csv", sep="")
    cc.name <- paste("ms02.", i, ".cc.csv", sep="")
    mf.name <- paste("ms02.", i, ".mf.csv", sep="")
    bp.mat <- matrix(as.numeric(unlist(read.table(bp.name, header=F, sep=","))), nrow=bp.dim, ncol=2)
    cc.mat <- matrix(as.numeric(unlist(read.table(cc.name, header=F, sep=","))), nrow=cc.dim, ncol=2)
    mf.mat <- matrix(as.numeric(unlist(read.table(mf.name, header=F, sep=","))), nrow=mf.dim, ncol=2)
    bp.perm <- rbind(bp.perm, c(bp.mat))
    cc.perm <- rbind(cc.perm, c(cc.mat))
    mf.perm <- rbind(mf.perm, c(mf.mat))
}

bp.mean <- colMeans(bp.perm)
bp.std <- colSds(bp.perm)
cc.mean <- colMeans(cc.perm)
cc.std <- colSds(cc.perm)
mf.mean <- colMeans(mf.perm)
mf.std <- colSds(mf.perm)

bp.zscore <- round((bp.obs - bp.mean)/bp.std, 3)
cc.zscore <- round((cc.obs - cc.mean)/cc.std, 3)
mf.zscore <- round((mf.obs - mf.mean)/mf.std, 3)

z.bp <- matrix(bp.zscore, nrow=bp.dim, ncol=2)
z.cc <- matrix(cc.zscore, nrow=cc.dim, ncol=2)
z.mf <- matrix(mf.zscore, nrow=mf.dim, ncol=2)

z.bp <- t(z.bp)
z.cc <- t(z.cc)
z.mf <- t(z.mf)

colnames(z.bp) <- bp.go.cat
rownames(z.bp) <- c("HMK", "Background")
colnames(z.cc) <- cc.go.cat
rownames(z.cc) <- c("HMK", "Background")
colnames(z.mf) <- mf.go.cat
rownames(z.mf) <- c("HMK", "Background")

bp.z.enrich <- z.bp[, which( z.bp[1,] >= 3 )]
bp.z.enrich <- bp.z.enrich[,order(-bp.z.enrich[1,])]
bp.enrich.list <- colnames(bp.z.enrich)
cc.z.enrich <- z.cc[, which( z.cc[1,] >= 3 )]
cc.z.enrich <- cc.z.enrich[,order(-cc.z.enrich[1,])]
cc.enrich.list <- colnames(cc.z.enrich)
mf.z.enrich <- z.mf[, which( z.mf[1,] >= 3 )]
mf.z.enrich <- mf.z.enrich[,order(-mf.z.enrich[1,])]
mf.enrich.list <- colnames(mf.z.enrich)

bp.z.suppress <- z.bp[, which( z.bp[1,] <= -3 )]
bp.z.suppress <- bp.z.suppress[,order(bp.z.suppress[1,])]
bp.suppress.list <- colnames(bp.z.suppress)
cc.z.suppress <- z.cc[, which( z.cc[1,] <= -3 )]
cc.z.suppress <- cc.z.suppress[,order(cc.z.suppress[1,])]
cc.suppress.list <- colnames(cc.z.suppress)
mf.z.suppress <- z.mf[, which( z.mf[1,] <= -3 )]
mf.z.suppress <- mf.z.suppress[,order(mf.z.suppress[1,])]
mf.suppress.list <- colnames(mf.z.suppress)

bp.enriched.terms <- c("GO.ID", "GO.Term", "HMK", "Background")
for (i in 1:length(bp.enrich.list)) {
  id <- as.character(bp.enrich.list[i])
  term <- Term(GOID(id))
  seri <- which(bp.go.cat %in% id)
  z.gene <- z.bp[1,seri]
  z.back <- z.bp[2,seri]
  bp.enriched.terms <- rbind(bp.enriched.terms, c(id, term, z.gene, z.back))
}
bp.enriched.terms
write.table(bp.enriched.terms, file="HMK.bp.enriched.txt", sep="\t", row.names = F, col.names = F, quote=T)

bp.suppressed.terms <- c("GO.ID", "GO.Term", "HMK", "Background")
for (i in 1:length(bp.suppress.list)) {
  id <- as.character(bp.suppress.list[i])
  term <- Term(GOID(id))
  seri <- which(bp.go.cat %in% id)
  z.gene <- z.bp[1,seri]
  z.back <- z.bp[2,seri]
  bp.suppressed.terms <- rbind(bp.suppressed.terms, c(id, term, z.gene, z.back))
}
bp.suppressed.terms
write.table(bp.suppressed.terms, file="HMK.bp.suppressed.txt", sep="\t", row.names = F, col.names = F, quote=T)

cc.enriched.terms <- c("GO.ID", "GO.Term", "HMK", "Background")
for (i in 1:length(cc.enrich.list)) {
  id <- as.character(cc.enrich.list[i])
  term <- Term(GOID(id))
  seri <- which(cc.go.cat %in% id)
  z.gene <- z.cc[1,seri]
  z.back <- z.cc[2,seri]
  cc.enriched.terms <- rbind(cc.enriched.terms, c(id, term, z.gene, z.back))
}
cc.enriched.terms
write.table(cc.enriched.terms, file="HMK.cc.enriched.txt", sep="\t", row.names = F, col.names = F, quote=T)

cc.suppressed.terms <- c("GO.ID", "GO.Term", "HMK", "Background")
for (i in 1:length(cc.suppress.list)) {
  id <- as.character(cc.suppress.list[i])
  term <- Term(GOID(id))
  seri <- which(cc.go.cat %in% id)
  z.gene <- z.cc[1,seri]
  z.back <- z.cc[2,seri]
  cc.suppressed.terms <- rbind(cc.suppressed.terms, c(id, term, z.gene, z.back))
}
cc.suppressed.terms
write.table(cc.suppressed.terms, file="HMK.cc.suppressed.txt", sep="\t", row.names = F, col.names = F, quote=T)

mf.enriched.terms <- c("GO.ID", "GO.Term", "HMK", "Background")
for (i in 1:length(mf.enrich.list)) {
  id <- as.character(mf.enrich.list[i])
  term <- Term(GOID(id))
  seri <- which(mf.go.cat %in% id)
  z.gene <- z.mf[1,seri]
  z.back <- z.mf[2,seri]
  mf.enriched.terms <- rbind(mf.enriched.terms, c(id, term, z.gene, z.back))
}
mf.enriched.terms
write.table(mf.enriched.terms, file="HMK.mf.enriched.txt", sep="\t", row.names = F, col.names = F, quote=T)

mf.suppressed.terms <- c("GO.ID", "GO.Term", "HMK", "Background")
for (i in 1:length(mf.suppress.list)) {
  id <- as.character(mf.suppress.list[i])
  term <- Term(GOID(id))
  seri <- which(mf.go.cat %in% id)
  z.gene <- z.mf[1,seri]
  z.back <- z.mf[2,seri]
  mf.suppressed.terms <- rbind(mf.suppressed.terms, c(id, term, z.gene, z.back))
}
mf.suppressed.terms
write.table(mf.suppressed.terms, file="HMK.mf.suppressed.txt", sep="\t", row.names = F, col.names = F, quote=T)

##Heatmaps
colors = c(seq(-2,2.4,length=10),seq(2.5,15.0,length=10))
my_palette <- colorRampPalette(c("blue2", "red2"))(n = 19)
##

bp.fig <- z.bp[, which( z.bp[1,] >= bp.z.enrich[1,10])]
bp.fig <- t(bp.fig[,order(-bp.fig[1,])])
png(filename = "bp.top10.png",width=3.5, height=6, res=1200, unit="in")
heatmap.2(bp.fig, col=my_palette, dendrogram='none', breaks=colors,
          sepcolor="lightgrey", sepwidth = c(0.05,0.01),
          colsep = 1:ncol(bp.z.enrich), rowsep = 1:nrow(bp.z.enrich),
          trace='none', Rowv=F, Colv=F, #revC=T, #labCol=NA,
          ylab="Biological Process Terms", xlab="",
          margins = c(3.2,5), key.title = NA, key.xlab=NA, key.ylab=NA,
          scale="none", symbreaks=F, symm=F, symkey = F,
          adjCol=c(0.5,0), adjRow=c(0.2,0), srtRow=45, srtCol=0,
          cexRow = 0.9, cexCol=0.9,
          cellnote = cbind(bp.fig[,1], bp.fig[,2]),
          main="")
dev.off()

mf.fig <- z.mf[, which( z.mf[1,] >= mf.z.enrich[1,10])]
mf.fig <- t(mf.fig[,order(-mf.fig[1,])])
png(filename = "mf.top10.png",width=3.5, height=6, res=1200, unit="in")
heatmap.2(mf.fig, col=my_palette, dendrogram='none', breaks=colors,
          sepcolor="lightgrey", sepwidth = c(0.05,0.01),
          colsep = 1:ncol(bp.z.enrich), rowsep = 1:nrow(bp.z.enrich),
          trace='none', Rowv=F, Colv=F, #revC=T, #labCol=NA,
          ylab="Molecular Function Terms", xlab="",
          margins = c(3.2,5), key.title = NA, key.xlab=NA, key.ylab=NA,
          scale="none", symbreaks=F, symm=F, symkey = F,
          adjCol=c(0.5,0), adjRow=c(0.2,0), srtRow=45, srtCol=0,
          cexRow = 0.9, cexCol=0.9,
          cellnote = cbind(mf.fig[,1], mf.fig[,2]),
          main="")
dev.off()

cc.fig <- z.cc[, which( z.cc[1,] >= cc.z.enrich[1,10])]
cc.fig <- t(cc.fig[,order(-cc.fig[1,])])
png(filename = "cc.top10.png",width=3.5, height=6, res=1200, unit="in")
heatmap.2(cc.fig, col=my_palette, dendrogram='none', breaks=colors,
          sepcolor="lightgrey", sepwidth = c(0.05,0.01),
          colsep = 1:ncol(bp.z.enrich), rowsep = 1:nrow(bp.z.enrich),
          trace='none', Rowv=F, Colv=F, #revC=T, #labCol=NA,
          ylab="Cellular Component Terms", xlab="",
          margins = c(3.2,5), key.title = NA, key.xlab=NA, key.ylab=NA,
          scale="none", symbreaks=F, symm=F, symkey = F,
          adjCol=c(0.5,0), adjRow=c(0.2,0), srtRow=45, srtCol=0,
          cexRow = 0.9, cexCol=0.9,
          cellnote = cbind(cc.fig[,1], cc.fig[,2]),
          main="")
dev.off()

