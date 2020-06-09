####
##This script lists DAVID and NetBAS result side-by-side
##p-values in DAVID are converted to z-scores, and scaled to match the maximum of NetBAS
library(pracma)
library(gplots)

filename <- read.csv("HMK.kegg.enriched.txt", header=T, sep="\t", stringsAsFactors=F)
pathway <- filename$Pathway[1:10]
netbas <- filename$NetBAS[1:10]
david <- filename$DAVID[1:10]

#convert DAVID p-values to z-scores
z.david <- c()
for (i in 1:10) {
    z.val <- abs(david[i] - sqrt(2)*erfcinv(2*david[i]))
    z.david <- rbind(z.david, z.val)
}

scale.factor <- max(netbas) / max(z.david[!is.na(z.david)])

z.david <- round(z.david * scale.factor,3)

zscores <- matrix(cbind(netbas, z.david), ncol=2)

colnames(zscores) <- c("NetPAS", "DAVID")
rownames(zscores) <- pathway

colors=c(seq(0,4.9,length=5), seq(5.1,max(netbas),length=5))
my_palette <- colorRampPalette(c("white", "red2"))(n=9)

p.david <- c("4.0E-106", "4.7E-83", "7.3E-68", "NA", "8.0E-52", "NA", "4.5E-73", "1.5E-15", "1.5E-6", "3.7E-3")
note <- cbind(netbas, p.david)

png(filename= "netbas-david.kegg.top10.png", width=5, height=6, res=1200, unit="in")
heatmap.2(zscores, col=my_palette, dendrogram='none', breaks=colors,
          colsep = 1:2, rowsep = 1:10, sepcolor="lightgrey", sepwidth=c(0.02,0.02),
          trace='none', Rowv=F, Colv=F, notecol="black",
          ylab="KEGG Pathways", xlab="",
          margins=c(1,12), key.title=NA, key.xlab=NA, key.ylab=NA,
          scale="none", symbreaks=F, symm=F, symkey=F,
          adjCol=c(0.5,0), adjRow=c(0.05,0), srtRow=45, srtCol=0,
          cexRow=0.9, cexCol=0.9, cellnote=note, main="")
dev.off()
