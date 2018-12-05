####
##This script lists DAVID and NetBAS result side-by-side
##p-values in DAVID are converted to z-scores, and scaled to match the maximum of NetBAS
rm(list=ls())
library(pracma)
library(gplots)

filename <- read.csv("HMK.cc.top10.txt", header=T, sep="\t", stringsAsFactors=F)
go.id <- filename$GO.ID
go.term <- filename$GO.Term
netbas <- filename$NetBAS
david <- filename$DAVID

#convert DAVID p-values to z-scores
z.david <- c()
for (i in 1:10) {
    z.val <- abs(david[i] - sqrt(2)*erfcinv(2*david[i]))
    z.david <- rbind(z.david, z.val)
}

scale.factor <- max(netbas) / max(z.david[!is.na(z.david)])

z.david <- round(z.david * scale.factor,3)

zscores <- matrix(cbind(netbas, z.david), ncol=2)

colnames(zscores) <- c("NetBAS", "DAVID")
rownames(zscores) <- go.id

colors=c(seq(0,4.9,length=5), seq(5.1,max(netbas),length=5))
my_palette <- colorRampPalette(c("white", "red2"))(n=9)

p.david <- c("1.1E-11", "NA", "1.7E-13", "2.4E-9", "5.9E-12", "NA", "3.8E-2", "NA", "NA", "1.1E-2")
note <- cbind(netbas, p.david)

png(filename= "netbas-david.cc.top10.png", width=3.5, height=6, res=1200, unit="in")
heatmap.2(zscores, col=my_palette, dendrogram='none', breaks=colors,
          colsep = 1:2, rowsep = 1:10, sepcolor="lightgrey", sepwidth=c(0.02,0.02),
          trace='none', Rowv=F, Colv=F,
          ylab="Cellular Component Terms", xlab="",
          margins=c(1,5.5), key.title=NA, key.xlab=NA, key.ylab=NA,
          scale="none", symbreaks=F, symm=F, symkey=F,
          adjCol=c(0.5,0), adjRow=c(0.05,0), srtRow=45, srtCol=0,
          cexRow=0.9, cexCol=0.9, cellnote=note, main="")
dev.off()
